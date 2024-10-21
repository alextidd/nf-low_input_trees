#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { preprocess                } from './modules/preprocess.nf'
include { reflag                    } from './modules/reflag.nf'
include { hairpin2                  } from './modules/hairpin2.nf'
include { post_filtering_and_pileup } from './modules/post_filtering_and_pileup.nf'
include { cgpVAF                    } from './modules/cgpVAF.nf'
include { sequoia                   } from './modules/sequoia.nf'

// print help message, return typical command line usage for the pipeline
if (params.help) {
  log.info paramsHelp("nextflow run low_input_trees --sample_sheet sample_sheet.csv --sequencing_type WGS --outdir out/")
  exit 0
}

// validate input parameters
if (params.validate_params) {
    validateParameters()
}

// Main workflow
workflow {

    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // preprocess
    preprocess()

    // reflag
    if ( params.reflag ) {
      reflag(preprocess.out.ch_input)
      hairpin2_ch = reflag.out
    }
    else {
      hairpin2_ch = preprocess.out.ch_input
    }

    // run hairpin
    hairpin2(hairpin2_ch)

    // run post-filtering and pileup
    post_filtering_and_pileup(hairpin2.out)

    // run cgpVAF
    cgpVAF(post_filtering_and_pileup.out, 
           preprocess.out.fasta, 
           preprocess.out.high_depth_bed,
           preprocess.out.cgpVAF_normal_bam)

    // run sequoia
    sequoia(cgpVAF.out,
            preprocess.out.fasta)
}

