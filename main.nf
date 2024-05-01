#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// PROCESSES
process hairpin_preselect {
    tag "${meta.sample_id}:${vcf_type}"
    label "normal"
    publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
      mode: "copy",
      pattern: "*.pre.vcf"
    publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
      mode: "symlink",
      pattern: "*.{annot.vcf.gz,sample.dupmarked.bam}"

    input:
    tuple val(meta), 
          val(vcf_type), path(vcf), 
          path(bam), path(bai), path(bas), path(met)

    output:
    tuple val(meta), 
          val(vcf_type), path(vcf), 
          path(bam), path(bai), path(bas), path(met),
          path("${meta.sample_id}.pre.vcf")

    script:
    """
    /bin/bash /code/runScriptPreselect.sh \
      --vcf-file ${vcf} \
      > ${meta.sample_id}.pre.vcf
    """
}

process hairpin_imitateANNOVAR {
    tag "${meta.sample_id}:${vcf_type}"
    label "normal"
    publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
      mode: "copy",
      pattern: "*.pre.annovar.txt"

    input:
    tuple val(meta), 
          val(vcf_type), path(vcf),
          path(bam), path(bai), path(bas), path(met),
          path(pre_vcf)

    output:
    tuple val(meta), 
          val(vcf_type), 
          path(bam), path(bai), path(bas), path(met),
          path(pre_vcf), path(vcf), 
          path("${meta.sample_id}.pre.annovar.txt")

    script:
    """
    /bin/bash /code/runScriptImitateANNOVAR.sh \
      --vcf-file ${pre_vcf} \
      > ${meta.sample_id}.pre.annovar.txt
    """
}

process hairpin_annotateBAMStatistics {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*.pre.annovar.annot.txt"

  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar)

  output:
  tuple val(meta), 
        val(vcf_type), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), path(vcf), 
        path(pre_annovar), 
        path("${meta.sample_id}.pre.annovar.annot.txt")

  script:
  """
  /bin/bash /code/runScriptAnnotate.sh \
    --annovarfile ${pre_annovar} \
    --bamfiles ${bam} \
    --threads ${task.cpus} \
    > ${meta.sample_id}.pre.annovar.annot.txt
  """
}

process hairpin_additionalBAMStatistics {
  tag "${meta.sample_id}:${vcf_type}"
  label "normal10gb"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/${meta.sample_id}", 
    mode: "copy",
    pattern: "*.pre.annovar.annot.full.txt"
  
  input:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot)
  path fasta
  path snp_database

  output:
  tuple val(meta), 
        val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path("${meta.sample_id}.pre.annovar.annot.full.txt")

  script:
  """
  /bin/bash /code/runScriptAdditional.sh \
    --annovarfile ${pre_annovar_annot} \
    --bamfile ${bam} \
    --threads $task.cpus \
    --reference ${fasta} \
    --snp-database ${snp_database} \
    > ${meta.sample_id}.pre.annovar.annot.full.txt
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
        path(pre_vcf), 
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
    --vcf-file ${pre_vcf} \
    --output-directory ./ \
    --prefix ${meta.sample_id} \
    --fragment-threshold ${params.fragment_threshold}
  """
}

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
        path("${meta.sample_id}_postfiltered.bed")

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
        val(sample_ids),
        val(vcf_type), path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(vcf_postfiltered)

  output:
  tuple val(meta),
        val(sample_ids),
        val(vcf_type), path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path("${meta.donor_id}_intervals.bed")

  script:
  """
  cut -f 1,2,4,5 *_postfiltered.vcf | sort | uniq \
  > ${meta.donor_id}_intervals.bed
  """
}

process cgpVAF_run {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"
  
  input: 
  tuple val(meta),
        val(sample_ids),
        val(vcf_type), path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals)
  tuple path(fasta), path(fai)
  tuple path(high_depth_bed), path(high_depth_tbi)
  tuple val(cgpVAF_normal_name), path(cgpVAF_normal_bam), path(cgpVAF_normal_bai)

  script:
  """
  module load cgpVAFcommand/2.5.0

  # get variant type
  [[ "$vcf_type" == "caveman" ]] && variant_type="snp"
  [[ "$vcf_type" == "pindel" ]] && variant_type="indel"

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
    --normal_name ${cgpVAF_normal_name} \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.join(' ')}
  """
}

process cgpVAFcommand {
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
  # module
  module purge
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

process get_mutation_filtering_parameters {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${run}", 
    mode: "copy"

  input:
  tuple val(meta),
        val(sample_ids),
        path(merged_caveman),
        path(merged_pindel)
        tuple val(run), val(extra_arguments)

  output:
  tuple val(meta),
        val(vcf_type),
        val(sample_ids)
  
  script:
  """
  get_mutation_filtering_parameters.R \
    --runid ${meta.experiment_id} \
    --snvfile ${merged_caveman} \
    --indelfile ${merged_pindel} \
    --output ./ ${extra_arguments}
  """
}

// SUBWORKFLOWS
workflow preprocess {
  main:
  // get fasta + fai
  Channel.fromPath(params.fasta, checkIfExists: true)
  | map { fasta -> [fasta, file(fasta + ".fai", checkIfExists: true)] } 
  | set { ch_fasta }

  // get bed + tbi
  Channel.fromPath(params.high_depth_bed, checkIfExists: true)
  | map { bed -> [bed, file(bed + ".tbi", checkIfExists: true)] }
  | set { ch_high_depth_bed }

  // get metadata + bam paths
  Channel.fromPath(params.sample_sheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      meta = [sample_id: row.sample_id, donor_id: row.donor_id, 
              project_id: row.project_id, experiment_id: row.experiment_id]
      [meta, 
      file(row.caveman_vcf, checkIfExists: true),
      file(row.pindel_vcf, checkIfExists: true),
      file(row.bam, checkIfExists: true),
      file(row.bam + ".bai", checkIfExists: true),
      file(row.bam + ".bas", checkIfExists: true),
      file(row.bam + ".met.gz", checkIfExists: true)]
  }
  | set { ch_samples }

  // get caveman vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta, "caveman", caveman_vcf, bam, bai, bas, met]
  }
  | set { ch_caveman }

  // get pindel vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta, "pindel", pindel_vcf, bam, bai, bas, met]
  }
  | set { ch_pindel }

  // concatenate channels together
  ch_caveman.concat(ch_pindel) 
  | set { ch_input } 

  emit:
  ch_input = ch_input
  ch_fasta = ch_fasta
  ch_high_depth_bed = ch_high_depth_bed
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

workflow post_filtering_and_pileup {
  take:
  ch_input

  main:
  // run post-filtering
  ch_input
  | post_filtering
  // group vcfs by donor, pileup
  | map { meta, vcf_type, vcf, vcf_postfiltered ->
          [meta.subMap('donor_id', 'project_id', 'experiment_id'), 
           meta.sample_id, vcf_type, vcfs, vcf_postfiltered] }
  | groupTuple(by: [0,1])
  | pileup

  emit:
  pileup.out
}

workflow cgpVAF {
  take:
  ch_input
  ch_fasta
  ch_high_depth_bed

  main:
  // run cpgVAF
  cgpVAFcommand(
    ch_input,
    ch_fasta,
    ch_high_depth_bed)

  emit:
  cgpVAFcommand.out
}

// Main workflow
workflow {

    // print help message, supply typical command line usage for the pipeline
    if (params.help) {
      log.info paramsHelp("nextflow run nf-chemo-trees --sample_sheet sample_sheet.csv")
      exit 0
    }

    // validate the input parameters
    // validateParameters()

    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // preprocess
    preprocess()

    // run hairpin
    hairpin(preprocess.out.ch_input)

    // run post_filtering_and_pileup
    post_filtering_and_pileup(hairpin.out)

    // run cgpVAF
    cgpVAF(post_filtering_and_pileup.out, 
           preprocess.out.ch_fasta, 
           preprocess.out.ch_high_depth_bed)

    // get mutation filtering parameters
    cgpVAF.out 
    | branch { meta, vcf_type, merged_muts ->
                caveman: vcf_type == "caveman"
                  return tuple(meta, merged_muts)
                pindel: vcf_type == "pindel"
                  return tuple(meta, merged_muts) }
    | set { ch_merged_muts }

    // run twice:
    // 1. not excluding any samples
    // 2. excluding samples with peak VAF < 0.4
    runs = [["include_all", ""],
            ["remove_below_0.4", "--mixed_remove --vaf_cutoff 0.4"]]
    ch_merged_muts.caveman
    | join(ch_merged_muts.pindel)
    | combine(runs)
    | view
}

