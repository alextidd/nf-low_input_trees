process pileup {
  tag "${meta.donor_id}:${meta.vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/", 
    mode: "copy",
    pattern: "*_intervals.bed"
  
  input:
  tuple val(meta),
        path(vcfs_filtered), path(tbis_filtered)

  output:
  tuple val(meta),
        path("${meta.donor_id}_intervals.bed")

  script:
  """
  (for file in *_filtered.vcf.gz ; do
    zgrep -v '^#' \$file | cut -f 1,2,4,5 ;
  done) | cat | sort | uniq \
  > ${meta.donor_id}_intervals.bed
  """
}

process cgpVAF_run {
  tag "${meta.donor_id}:${meta.vcf_type}"
  label "normal" // label "week50gb"
  maxRetries 3

  input: 
  tuple val(meta),
        val(sample_ids),
        path(vcfs), path(tbis),
        path(bams), path(bais), path(bass), path(mets),
        path(intervals_bed)
  tuple path(fasta), path(fai)
  tuple path(cgpvaf_high_depth_bed), path(high_depth_tbi)
  tuple path(cgpvaf_normal_bam),
        path(cgpVAF_normal_bas),
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),   
        path("${meta.donor_id}_${meta.vcf_type}_vaf.tsv")

  script:
  def variant_type = (meta.vcf_type == "caveman") ? "snp" : (meta.vcf_type == "pindel") ? "indel" : ""
  """
  # modules
  module load cgpVAFcommand/2.5.0

  # run cgpVAF
  cgpVaf.pl \\
    --inputDir ./ \\
    --outDir ./ \\
    --variant_type ${variant_type} \\
    --genome ${fasta} \\
    --cgpvaf_high_depth_bed ${cgpvaf_high_depth_bed} \\
    --bed_only 1 \\
    --bedIntervals ${intervals_bed} \\
    --base_quality 25 \\
    --map_quality 30 \\
    --tumour_name ${sample_ids.join(' ')} \\
    --tumour_bam ${bams.join(' ')} \\
    --normal_name "normal" \\
    --normal_bam ${cgpvaf_normal_bam} \\
    --vcf ${vcfs.join(' ')}

  # rename output and remove vcf header
  cgpVAF_out=(`ls *_vaf.tsv`)
  n_files=\${#cgpVAF_out[@]}
  if (( n_files > 1 )) ; then
    echo "Error: More than one file matches pattern *_vaf.tsv"
    exit 1
  elif (( n_files == 1 )) ; then
    grep -v '^##' \${cgpVAF_out[0]} \
    > ${meta.donor_id}_${meta.vcf_type}_vaf.tsv
  else
    echo "Error: No files match the pattern *_vaf.tsv"
    exit 1
  fi 
  """
}

process cgpVAF_merge {
  tag "${meta.donor_id}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/", mode: 'copy' 

  input:
  tuple val(meta),
        path(vafs, stageAs: "vaf??.tsv")

  output:
  tuple val(meta),
        path("${meta.donor_id}_${meta.vcf_type}_vaf.tsv")

  script:
  """
  # get variant info cols: Chrom Pos Ref Alt normal_MTR normal_DEP
  var_cols=\$(head -1 ${vafs[0]} \\
    | tr '\\t' '\\n' \\
    | grep -n 'Chrom\\|Pos\\|Ref\\|Alt\\|normal_MTR\\|normal_DEP' \\
    | cut -d: -f1 | tr '\\n' ',' | sed 's/,\$//')
  cut -f\$var_cols ${vafs[0]} \\
  > tmp.1.cut

  # for all files, get *_MTR and *_DEP columns
  for file in ${vafs.join(' ')} ; do
    cols=\$(head -1 \$file | tr '\\t' '\\n' \\
      | grep -n '_\\(MTR\\|DEP\\)\$' \\
      | grep -v '^.*normal_\\(DEP\\|MTR\\)\$' \\
      | cut -d: -f1 | tr '\\n' ',' | sed 's/,\$//')
    cut -f\$cols \$file > tmp.\$file.cut
  done

  # combine all extracted columns
  paste tmp.*.cut > ${meta.donor_id}_${meta.vcf_type}_vaf.tsv
  """
}

workflow cgpVAF {
  take:
  ch_input
  fasta
  cgpvaf_high_depth_bed
  cgpvaf_normal_bam

  main:

  // group vcfs by donor, pileup
  ch_input
  | map { meta, vcf_filtered, tbi_filtered, bam, bai, bas, met  ->
          [meta.subMap("donor_id", "vcf_type"), vcf_filtered, tbi_filtered] }
  | groupTuple
  | pileup

  // batch samples per donor
  ch_input
  | map { meta, vcf_filtered, tbi_filtered, bam, bai, bas, met  ->
          [meta.subMap("donor_id", "vcf_type"), 
           meta.sample_id, vcf_filtered, tbi_filtered, bam, bai, bas, met] }
  | groupTuple(size: params.cgpvaf_sample_batch_size, remainder: true)
  | combine(pileup.out, by: 0)
  | set { ch_batched_samples }

  // run batched cgpvaf and then merge
  // TODO: move bgzip and tabix of filtered vcfs outside cgpVAF_run
  cgpVAF_run(ch_batched_samples, fasta, cgpvaf_high_depth_bed, cgpvaf_normal_bam)
  | groupTuple
  | cgpVAF_merge
  | map { meta, cgpVAF_out -> 
          [meta.subMap("donor_id"), cgpVAF_out]}
  | groupTuple
  | set { ch_cgpVAF_out }

  emit:
  ch_cgpVAF_out
}