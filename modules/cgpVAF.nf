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
        path("${meta.donor_id}_merged_${vcf_type}.tsv")

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

  # create full matrix, concat donor results
  vaf_files=\$(ls *_vaf.tsv)

  # get positional columns
  cat \${vaf_files[0]} | grep -v '^##' | 
  cut -f3,4,5,6,24,26 \
  > ${meta.donor_id}_merged_${vcf_type}.positions.tsv.tmp

  # get variant columns
  for file in \${vaf_files[@]} ; do
    cat \$file | 
    grep -v '^##' |
    cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 \
    > \${file}.tmp
  done

  # paste together
  paste ${meta.donor_id}_merged_${vcf_type}.positions.tsv.tmp *_vaf.tsv.tmp \
  > ${meta.donor_id}_merged_${vcf_type}.tsv
  """
}

process cgpVAF_command {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"

  input:
  tuple val(meta),
        val(vcf_type), path(vcfs), 
        val(sample_ids),
        path(bed)
  tuple path(fasta), path(fai)
  tuple path(high_depth_bed), path(high_depth_tbi)

  output:
  tuple val(meta),
        val(vcf_type),
        path("${meta.donor_id}_merged_${vcf_type}.tsv")

  script:
  """
  module load cgpVAFcommand/2.5.0

  # get input code for vcf type
  [[ "$vcf_type" == "pindel" ]] && in="1"
  [[ "$vcf_type" == "caveman" ]] && in="3"

  # generate cgpVAF command
  echo \$in |
  createVafCmd.pl \
    -pid ${meta.project_id} \
    --outdir ./ \
    --genome ${fasta} \
    --high_depth_bed ${high_depth_bed} \
    -bq 25 -mq 30 -bo 1 \
    --bedIntervals ${bed} \
    --sample_names ${sample_ids.join(' ')}

  # run cgpVAF
  bash cgpVafChr.cmd 

  # concat chromosomal results
  bash cgpVafConcat.cmd 

  # create full matrix, concat donor results
  vaf_files=\$(ls output/*/*/*_vaf.tsv)

  # get positional columns
  cat \${vaf_files[0]} | grep -v '^##' | 
  cut -f3,4,5,6,24,26 \
  > ${meta.donor_id}_merged_${vcf_type}.tsv.tmp

  # get variant columns
  for file in \${vaf_files[@]} ; do
    cat \$file | 
    grep -v '^##' |
    cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 \
    > \${file}.tmp
  done

  # paste together
  paste ${meta.donor_id}_merged_${vcf_type}.tsv.tmp output/*/*/*_vaf.tsv.tmp \
  > ${meta.donor_id}_merged_${vcf_type}.tsv
  """
}

workflow cgpVAF {
  take:
  ch_input
  ch_fasta
  ch_high_depth_bed
  ch_cgpVAF_normal_bam

  main:
  // run cpgVAF
  //cgpVAF_command(
  //  ch_input,
  //  ch_fasta,
  //  ch_high_depth_bed)

  cgpVAF_run(
    ch_input,
    ch_fasta,
    ch_high_depth_bed,
    ch_cgpVAF_normal_bam
  )

  emit:
  cgpVAF_run.out
}