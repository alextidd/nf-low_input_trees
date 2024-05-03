#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'
include { preprocess } from './modules/preprocess.nf'
include { hairpin } from './modules/hairpin.nf'
include { post_filtering_and_pileup } from './modules/post_filtering_and_pileup.nf'
include { cgpVAF } from './modules/cgpVAF.nf'
include { sequoia } from './modules/sequoia.nf'

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
    hairpin(preprocess.out.ch_input,
            preprocess.out.genome_names)

    // run post-filtering and pileup
    post_filtering_and_pileup(hairpin.out)

    // run cgpVAF
    cgpVAF(post_filtering_and_pileup.out, 
           preprocess.out.ch_fasta, 
           preprocess.out.ch_high_depth_bed,
           preprocess.out.ch_cgpVAF_normal_bam)

    // run sequoia
    sequoia(cgpVAF.out,
            preprocess.out.ch_fasta)
    
}

