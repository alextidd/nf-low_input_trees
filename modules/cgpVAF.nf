process pileup {
  tag "${meta.donor_id}:${meta.vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/", 
    mode: "copy",
    pattern: "*_intervals.bed"
  
  input:
  tuple val(meta),
        path(vcfs_postfiltered)

  output:
  tuple val(meta),
        path("${meta.donor_id}_intervals.bed")

  script:
  """
  (for file in *_postfiltered.vcf ; do
    grep -v '^#' \$file | cut -f 1,2,4,5 ;
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
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(intervals_bed)
  tuple path(fasta), path(fai)
  tuple path(high_depth_bed), path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam),
        path(cgpVAF_normal_bas),
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),   
        path("${meta.donor_id}_${meta.vcf_type}_vaf.tsv")

  script:
  def variant_type = (meta.vcf_type == "caveman") ? "snp" : (meta.vcf_type == "pindel") ? "indel" : ""
  """
  module load cgpVAFcommand/2.5.0

  # gzip the vcfs
  for file in *.vcf ; do
    gzip -f \${file}
  done

  # run cgpVAF
  cgpVaf.pl \
    --inputDir ./ \
    --outDir ./ \
    --variant_type ${variant_type} \
    --genome ${fasta} \
    --high_depth_bed ${high_depth_bed} \
    --bed_only 1 \
    --bedIntervals ${intervals_bed} \
    --base_quality 25 \
    --map_quality 30 \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.collect { it + '.gz' }.join(' ')}

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

workflow cgpVAF {
  take:
  ch_input
  fasta
  high_depth_bed
  cgpVAF_normal_bam

  main:

  // group vcfs by donor, pileup
  ch_input
  | map { meta, vcf_postfiltered, bam, bai, bas, met  ->
          [meta.subMap("donor_id", "vcf_type"), vcf_postfiltered] }
  | groupTuple
  | pileup

  // batch samples per donor
  ch_input
  | map { meta, vcf_postfiltered, bam, bai, bas, met  ->
          [meta.subMap("donor_id", "vcf_type"), 
           meta.sample_id, vcf_postfiltered, bam, bai, bas, met] }
  | groupTuple(size: params.cgpVAF_sample_batch_size, remainder: true)
  | combine(pileup.out, by: 0)
  | set { ch_batched_samples }

  // run cgpVAF, group outputs by donor
  cgpVAF_run(ch_batched_samples, fasta, high_depth_bed, cgpVAF_normal_bam)
  | map { meta, cgpVAF_out -> 
          [meta.subMap("donor_id"), cgpVAF_out]}
  | groupTuple
  | set { ch_cgpVAF_out }

  ch_cgpVAF_out.view()

  emit:
  ch_cgpVAF_out
}