process filtering_run {
  tag "${meta.sample_id}:${meta.vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*_postfiltered.vcf"

  input:
  tuple val(meta),
        path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta),
        path("${meta.sample_id}_postfiltered.vcf"),
        path(bam), path(bai), path(bas), path(met)

  script:
  if (meta.vcf_type == "caveman") {
    """
    # apply filters as described in https://confluence.sanger.ac.uk/display/CAS/hairpin
    # FILTER = PASS
    # INFO/CLPM=0.00 (A soft flag median number of soft clipped bases in variant supporting reads)
    # INFO/ASRD>=0.87 (A soft flag median (read length adjusted) alignment score of reads showing the variant allele)

    module load bcftools-1.9/python-3.11.6
    bcftools filter \\
      -i 'FILTER="PASS" && INFO/CLPM="0.00" && INFO/ASRD>=0.87' \\
      ${vcf} \\
    > ${meta.sample_id}_postfiltered.vcf
    """
  } else if (meta.vcf_type == "pindel") {
    """
    # apply pass filter (FILTER = PASS)

    module load bcftools-1.9/python-3.11.6
    bcftools filter \\
      -i 'FILTER="PASS"' \\
      ${vcf} \\
    > ${meta.sample_id}_postfiltered.vcf
    """
  }
}

workflow filtering {
  take:
  ch_input

  main:
  // run post-filtering
  ch_input
  | filtering_run

  emit:
  filtering_run.out
}