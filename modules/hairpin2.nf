// filter VCF on FILTER=PASS and CLPM=0 and ASMD > 140
process hairpin2 {
  tag "${meta.sample_id}:${meta.vcf_type}"
  label "normal"

  input:
  tuple val(meta), 
        path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}_hairpin.vcf")

  script:
  """
  zcat ${vcf} > ${meta.sample_id}_reflagged.vcf
  hairpin2 \\
    --vcf-in ${meta.sample_id}_reflagged.vcf \\
    --vcf-out ${meta.sample_id}_hairpin.vcf \\
    --alignments ${bam} \\
    --format b \\
    --name-mapping TUMOUR:${meta.sample_id}
  """
}


process sync_pindels {
  tag "${meta.sample_id}:${meta.vcf_type}"
  label "normal"
  
  input:
  tuple val(meta), 
        path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}_hairpin.vcf")

  script:
  """
  zcat $vcf > ${meta.sample_id}_hairpin.vcf
  """
}

workflow hairpin {
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
  | hairpin2

  // get pindel channel into same format as hairpin2 output
  // (meta, vcf, bam, bai, bas, met, vcf_hairpin2)
  // (vcf_hairpin2 here is the same file as vcf, but unzipped and renamed)
  ch_branched.pindel
  | sync_pindels

  sync_pindels.out.concat(hairpin2.out) 
  | set { ch_hairpin2 } 

  emit:
  ch_hairpin2
}


