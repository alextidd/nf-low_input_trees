// filter VCF on FILTER=PASS and CLPM=0 and ASMD > 140
process hairpin2_run {
  tag "${meta.sample_id}:${meta.vcf_type}"
  label "normal"
  errorStrategy 'retry'

  input:
  tuple val(meta), 
        path(vcf),
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        path("${meta.sample_id}_hairpin2.vcf"),
        path(bam), path(bai), path(bas), path(met)
  
  script:
  """
  # modules
  module load hairpin2/hairpin2-1.0.0
  module load samtools-1.19/python-3.12.0

  # index
  tabix -f -p vcf ${vcf}
  
  # run hairpin2
  hairpin2 \\
    --vcf-in ${vcf} \\
    --vcf-out ${meta.sample_id}_hairpin2.vcf \\
    --alignments ${bam} \\
    --format b \\
    --name-mapping TUMOUR:${meta.sample_id}
  """
}

workflow hairpin2 {
  take: 
  ch_input

  main:

  // branch caveman and pindel outputs
  ch_input
  | branch {
      meta, vcf, bam, bai, bas, met ->
      caveman: meta.vcf_type == "caveman"
        return tuple(meta, vcf, bam, bai, bas, met)
      pindel: meta.vcf_type == "pindel"
        return tuple(meta, vcf, bam, bai, bas, met)
  } | set { ch_branched }

  // run hairpin on caveman
  ch_branched.caveman 
  | hairpin2_run

  // concat channels
  ch_branched.pindel.concat(hairpin2_run.out) 
  | set { ch_hairpin2 } 

  emit:
  ch_hairpin2
}


