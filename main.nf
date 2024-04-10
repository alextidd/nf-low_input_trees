#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

process preselect {
    tag "${meta.donor_id}:${meta.sample_id}"
    label "normal10gb"

    input:
    tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met)

    output:
    tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path("${meta.sample_id}.pre.vcf")

    script:
    """
    module load singularity
    singularity run \
    --bind /nfs,/lustre \
    --app preselect ${params.sif} \
    --vcf-file ${vcf} \
    > ${meta.sample_id}.pre.vcf
    """
}

process imitateANNOVAR {
    tag "${meta.donor_id}:${meta.sample_id}"
    label "normal10gb"

    input:
    tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf)

    output:
    tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path("${meta.sample_id}.pre.annovar.txt")

    script:
    """
    module load singularity
    singularity run \
    --bind /nfs,/lustre \
    --app imitateANNOVAR ${params.sif} \
    --vcf-file ${pre_vcf} \
    > ${meta.sample_id}.pre.annovar.txt
    """
}

process annotateBAMStatistics {
  tag "${meta.donor_id}:${meta.sample_id}"
  label "normal10gb"

  input:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar)

  output:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path("${meta.sample_id}.pre.annovar.annot.txt")

  script:
  """
  module load singularity
  singularity run \
    --bind /nfs,/lustre \
    --app annotateBAMStatistics ${params.sif} \
    --annovarfile ${pre_annovar} \
    --bamfiles ${bam} \
    --threads ${task.cpus} \
    > ${meta.sample_id}.pre.annovar.annot.txt
  """
}

process additionalBAMStatistics {
  tag "${meta.donor_id}:${meta.sample_id}"
  label "normal10gb"

  input:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot)

  output:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path("${meta.sample_id}.pre.annovar.annot.full.txt")

  script:
  """
  module load singularity
  singularity run \
    --bind /nfs,/lustre \
    --app additionalBAMStatistics ${params.sif} \
    --annovarfile ${pre_annovar_annot} \
    --bamfile ${bam} \
    --threads $task.cpus \
    --reference ${params.fasta} \
    --snp-database ${params.snp_database} \
    > ${meta.sample_id}.pre.annovar.annot.full.txt
  """
}

process filtering {
  tag "${meta.donor_id}:${meta.sample_id}"
  label "normal10gb"

  input:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path(pre_annovar_annot_full)

  output:
  tuple val(meta), val(vcf_type), path(vcf), 
        path(bam), path(bai), path(bas), path(met),
        path(pre_vcf), 
        path(pre_annovar), 
        path(pre_annovar_annot), 
        path(pre_annovar_annot_full), 
        path("${meta.sample_id}.pre.annovar.annot.full.txt")

  script:
  """
  module load singularity
  singularity run \
    --bind /nfs,/lustre \
    --app filtering ${params.sif} \
    --annotated-file ${pre_annovar_annot_full} \
    --vcf-file ${pre_vcf} \
    --output-dir ./ \
    --prefix ${meta.sample_id} \
    --fragment-threshold ${params.fragment_threshold}
  """
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

    // create a new channel of metadata from a sample sheet
    //ch_input = Channel.fromSamplesheet("sample_sheet")
    //ch_input.view()

    // get metadata + bam paths
    Channel.fromPath(params.sample_sheet, checkIfExists: true)
    | splitCsv(header: true)
    | map { row ->
        meta = [sample_id: row.sample_id, donor_id: row.donor_id]
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
    ch_samples \
    | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
            [meta, "pindel", pindel_vcf, bam, bai, bas, met]
    }
    | set { ch_pindel }

    // concatenate channels together
    ch_caveman.concat(ch_pindel) 
    | set { ch_input } 

    // run
    ch_input 
    | preselect
    | imitateANNOVAR
    | annotateBAMStatistics
    | additionalBAMStatistics
    | filtering

}

