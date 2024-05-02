
process post_filtering {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  errorStrategy 'ignore'

  input:
  tuple val(meta),
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(passed_vcf),
        path(filtered_vcf)

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
      ${passed_vcf} |
    grep -v '^#' \
    > ${meta.sample_id}_postfiltered.vcf
    """
  } else if (vcf_type == "pindel") {
    """
    # apply pass filter (FILTER = PASS)

    module load bcftools-1.9/python-3.11.6
    bcftools filter \
      -i 'FILTER="PASS"' \
      ${passed_vcf} |
    grep -v '^#' \
    > ${meta.sample_id}_postfiltered.vcf
    """
  }
}

process pileup {
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
        path(vcf_postfiltered)

  output:
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path("${meta.donor_id}_intervals.bed")

  script:
  """
  cut -f 1,2,4,5 *_postfiltered.vcf | sort | uniq \
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
          [meta.subMap('donor_id', 'project_id', 'experiment_id'), 
           vcf_type, meta.sample_id, vcf, bam, bai, bas, met, vcf_postfiltered] }
  | groupTuple(by: [0,1])
  | pileup

  emit:
  pileup.out
}
