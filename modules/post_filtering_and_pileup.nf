process post_filtering {
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

process pileup {
  tag "${meta.donor_id}:${meta.vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/", 
    mode: "copy",
    pattern: "*_intervals.bed"
  
  input:
  tuple val(meta),
        val(sample_ids),
        path(vcfs_postfiltered), 
        path(bams), path(bais), path(bass), path(mets)

  output:
  tuple val(meta),
        val(sample_ids),
        path(vcfs_postfiltered), 
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
  | map { meta, vcf_postfiltered, bam, bai, bas, met  ->
          [meta.subMap("donor_id", "vcf_type"), 
           meta.sample_id, vcf_postfiltered, bam, bai, bas, met] }
  | groupTuple
  | pileup

  emit:
  pileup.out
}