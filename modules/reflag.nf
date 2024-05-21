process reflag_pindel {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  script:
  """
  
  """
}

process reflag_caveman {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  script:
  """
  
  """
}

workflow hairpin {
  take: 
  ch_caveman
  ch_pindel

  main:
  // reflag caveman
  ch_caveman 
  | reflag_caveman

  // reflag pindel
  ch_pindel
  | reflag_pindel

  emit:
  reflag_caveman.out
  reflag_pindel.out
}


