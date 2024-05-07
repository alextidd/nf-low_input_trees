// filter VCF on FILTER=PASS and CLPM=0 and ASMD > 140
process hairpin_preselect {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "symlink",
    pattern: "*.{annot.vcf.gz,sample.dupmarked.bam}"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*.preselected.vcf"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}.preselected.vcf")

  script:
  """
  # disable ASMD filter (will filter on ASRD in the next step)
  /bin/bash /code/runScriptPreselect.sh \
    --vcf-file ${vcf} \
    > ${meta.sample_id}.preselected.vcf
  """
}

// convert the VCF file to ANNOVAR format (Chr,Start,End,Ref,Alt)
process hairpin_imitateANNOVAR {
    tag "${meta.sample_id}:${vcf_type}"
    label "normal"
    publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
      mode: "copy",
      pattern: "*.preselected.annovar.txt"

    input:
    tuple val(meta), 
          val(vcf_type), path(vcf),
          path(bam), path(bai), path(bas), path(met),
          path(preselected_vcf)

    output:
    tuple val(meta), 
          val(vcf_type), path(vcf),
          path(bam), path(bai), path(bas), path(met),
          path(preselected_vcf), 
          path("${meta.sample_id}.preselected.annovar.txt")

    script:
    """
    /bin/bash /code/runScriptImitateANNOVAR.sh \
      --vcf-file ${preselected_vcf} \
      > ${meta.sample_id}.preselected.annovar.txt
    """
}

// 
process hairpin_annotateBAMStatistics {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*.preselected.annovar.annot.txt"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(preselected_vcf), 
        path(pre_annovar)

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(preselected_vcf), 
        path(pre_annovar), 
        path("${meta.sample_id}.preselected.annovar.annot.txt")

  script:
  """
  /bin/bash /code/runScriptAnnotate.sh \
    --annovarfile ${pre_annovar} \
    --bamfiles ${bam} \
    --threads ${task.cpus} \
    > ${meta.sample_id}.preselected.annovar.annot.txt
  """
}

process hairpin_additionalBAMStatistics {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal10gb"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*.preselected.annovar.annot.full.txt"
  
  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(preselected_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot)
  path fasta
  path snp_database

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(preselected_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path("${meta.sample_id}.preselected.annovar.annot.full.txt")

  script:
  """
  /bin/bash /code/runScriptAdditional.sh \
    --annovarfile ${pre_annovar_annot} \
    --bamfile ${bam} \
    --threads $task.cpus \
    --reference ${fasta} \
    --snp-database ${snp_database} \
    > ${meta.sample_id}.preselected.annovar.annot.full.txt
  """
}

process hairpin_filtering {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*{passed,filtered}.vcf"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(preselected_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path(pre_annovar_annot_full)

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}_passed.vcf"),
        path("${meta.sample_id}_filtered.vcf")

  script:
  """
  /bin/bash /code/runScriptFiltering.sh \
    --annotated-file ${pre_annovar_annot_full} \
    --vcf-file ${preselected_vcf} \
    --output-directory ./ \
    --prefix ${meta.sample_id} \
    --fragment-threshold ${params.fragment_threshold}
  """
}

workflow hairpin {
  take: 
  ch_input

  main:
  // run
  ch_input 
  | hairpin_preselect
  | hairpin_imitateANNOVAR
  | hairpin_annotateBAMStatistics

  // add fasta and snp database to input
  hairpin_additionalBAMStatistics (
    hairpin_annotateBAMStatistics.out,
    params.fasta,
    params.snp_database)
  | hairpin_filtering

  emit:
  hairpin_filtering.out
}


