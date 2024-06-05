process post_filtering {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*_postfiltered.vcf"

  input:
  tuple val(meta),
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(vcf_passed, stageAs: "passed.vcf")

  output:
  tuple val(meta),
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}_postfiltered.vcf")

  script:
  if (vcf_type == "caveman") {
    """
    # apply filters as described in https://confluence.sanger.ac.uk/display/CAS/hairpin
    # FILTER = PASS
    # INFO/CLPM=0.00 (A soft flag median number of soft clipped bases in variant supporting reads)
    # INFO/ASRD>=0.87 (A soft flag median (read length adjusted) alignment score of reads showing the variant allele)

    module load bcftools-1.9/python-3.11.6
    bcftools filter \
      -i 'FILTER="PASS" && INFO/CLPM="0.00" && INFO/ASRD>=0.87' \
      ${vcf_passed} \
    > ${meta.sample_id}_postfiltered.vcf
    """
  } else if (vcf_type == "pindel") {
    """
    # apply pass filter (FILTER = PASS)

    module load bcftools-1.9/python-3.11.6
    bcftools filter \
      -i 'FILTER="PASS"' \
      ${vcf_passed} \
    > ${meta.sample_id}_postfiltered.vcf
    """
  }
}

process pileup {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy",
    pattern: "*_intervals.bed"
  
  input:
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(vcfs_postfiltered)

  output:
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path("${meta.donor_id}_intervals.bed")

  script:
  """
  (for file in *_postfiltered.vcf ; do
    grep -v '^#' \$file | cut -f 1,2,4,5 ;
  done) | cat | sort | uniq \
  > ${meta.donor_id}_intervals.bed
  """
}

workflow post_filtering_and_pileup {
  take:
  ch_input

  main:
  // run post-filtering
  ch_input
  | post_filtering
  // group vcfs by donor, pileup
  | map { meta, vcf_type, vcf, bam, bai, bas, met, vcf_postfiltered ->
          [meta.subMap('donor_id'), 
           vcf_type, meta.sample_id, vcf, bam, bai, bas, met, vcf_postfiltered] }
  | groupTuple(by: [0,1])
  | pileup

  // ch_input
  // | post_filtering
  // | map { meta, vcf_type, vcf, bam, bai, bas, met, vcf_postfiltered ->
  //         [meta.subMap('donor_id'), 
  //          vcf_type, meta.sample_id, vcf, bam, bai, bas, met, vcf_postfiltered] }
  // | set { ch_to_batch }

  // // group vcfs by donor, pileup intervals
  // ch_to_batch
  // | groupTuple(by: [0,1])
  // | pileup

  // // batch vcfs / bams by 10 within donors, combine with bed intervals
  // ch_to_batch
  // | flatMap { meta, vcf_type, vcf, bam, bai, bas, met, vcf_postfiltered -> 
  //             vcf.collate(50).collect { chunk -> [id, chunk] } }

  emit:
  pileup.out
}
