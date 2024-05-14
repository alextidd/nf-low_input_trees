process cgpVAF_run {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"
  
  input: 
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals)
  tuple path(fasta), 
        path(fai)
  tuple path(high_depth_bed), 
        path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),
        path("${meta.donor_id}_${vcf_type}_vaf.tsv")

  script:
  """
  module load cgpVAFcommand/2.5.0

  # get variant type
  [[ "$vcf_type" == "caveman" ]] && variant_type="snp"
  [[ "$vcf_type" == "pindel" ]] && variant_type="indel"

  # run cgpVAF
  cgpVaf.pl \
    --inputDir ./ \
    --outDir ./ \
    --variant_type \$variant_type \
    --genome ${fasta} \
    --high_depth_bed ${params.high_depth_bed} \
    --bed_only 1 \
    --bedIntervals ${bed_intervals} \
    --base_quality 25 \
    --map_quality 30 \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.join(' ')}
  
  # rename output and remove vcf header
  cgpVAF_out=(`ls *_vaf.tsv`)
  n_files=\${#cgpVAF_out[@]}
  if (( n_files > 1 )) ; then
    echo "Error: More than one file matches pattern *_vaf.tsv"
    exit 1
  elif (( n_files == 1 )) ; then
    grep -v '^##' \${cgpVAF_out[0]} \
    > ${meta.donor_id}_${vcf_type}_vaf.tsv
  else
    echo "Error: No files match the pattern *_vaf.tsv"
    exit 1
  fi 
  """
}

workflow cgpVAF {
  take:
  ch_input
  ch_fasta
  ch_high_depth_bed
  ch_cgpVAF_normal_bam

  main:
  cgpVAF_run(
    ch_input,
    ch_fasta,
    ch_high_depth_bed,
    ch_cgpVAF_normal_bam)
  | groupTuple
  | set { cgpVAF_out }

  emit:
  cgpVAF_out
}