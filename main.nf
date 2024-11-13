#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; samplesheetToList } from 'plugin/nf-schema'
include { preprocess                } from './modules/preprocess.nf'
include { reflag                    } from './modules/reflag.nf'
include { hairpin2                  } from './modules/hairpin2.nf'
include { filtering                 } from './modules/filtering.nf'
include { cgpVAF                    } from './modules/cgpVAF.nf'
include { sequoia                   } from './modules/sequoia.nf'

// Main workflow
workflow {

    // print help message, return typical command line usage for the pipeline
    if (params.help) {
      log.info paramsHelp("nextflow run low_input_trees --sample_sheet sample_sheet.csv --sequencing_type WGS --outdir out/")
      exit 0
    }

    // validate input params
    if (params.validate_params) {
        validateParameters()
    }

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

    // run hairpin2
    hairpin2(hairpin2_ch)

    // run post-filtering and pileup
    filtering(hairpin2.out)

    // run cgpVAF
    cgpVAF(filtering.out,
           preprocess.out.fasta,
           preprocess.out.cgpvaf_high_depth_bed,
           preprocess.out.cgpvaf_normal_bam)

    // run sequoia
    sequoia(cgpVAF.out,
            preprocess.out.fasta)
}

