#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

//export NXF_PLUGINS_DIR=~/.nextflow/plugins/

// all of the default parameters are being set in `nextflow.config`

// import sub-workflows
include { validate_manifest } from './modules/manifest'
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// Function which prints help message text
def helpMessage() {
    log.info"""
Usage:

nextflow run nf-chemotrees <ARGUMENTS>

Required Arguments:

  Input Data:
  --sample_sheet        Single file with the location of all input data. Must be formatted
                        as a CSV containing columns: sample_id,donor_id,bam_file,...

  Reference Data:
  --genome_fasta        Reference genome to use for alignment, in FASTA format

  Output Location:
  --out_dir             Output directory

Optional Arguments:
  --...
    """.stripIndent()
}


// Main workflow
workflow {

    // Show help message if the user specifies the --help flag at runtime
    // or if any required params are not provided
    if ( params.help || params.output_folder == false || params.genome_fasta == false ){
        // Invoke the function above which prints the help message
        helpMessage()
        // Exit out and do not run anything else
        exit 1
    }

    // validate the input parameters
    validateParameters()

    // print summary of supplied parameters
    log.info paramsSummaryLog(workflow)

    // create a new channel of metadata from a sample sheet
    ch_input = Channel.fromSamplesheet(params.sample_sheet)
    ch_input.view()
    

}