#!/usr/bin/env nextflow

// using DSL-2
nextflow.enable.dsl=2

// all of the default parameters are being set in `nextflow.config`

// import functions / modules / subworkflows / workflows
include { validateParameters; paramsHelp; paramsSummaryLog; fromSamplesheet } from 'plugin/nf-validation'

// run hairpin
process HAIRPIN {
  tag "${meta.donor_id}:${meta.sample_id}"
  label "normal10gb"

  input:
  tuple val(meta), path(caveman_vcf), path(pindel_vcf), path(bam), path(bai), path(bas), path(met)
  path fasta

  script:
  """
  # module load hairpin
  # hairpin \
  #   -v ${caveman_vcf} \
  #   -b ${bam}
  
  # run singularity images
  singularity pull shub://MathijsSanders/SangerLCMFilteringSingularity
  singularity run --bind /nfs,/lustre --app preselect SangerLCMFilteringSingularity_latest.sif -v ${caveman_vcf} > filtered.vcf
  singularity run --bind /nfs,/lustre --app imitateANNOVAR SangerLCMFilteringSingularity_latest.sif -v filtered.vcf > filtered.annovar
  singularity run --bind /nfs,/lustre --app annotateBAMStatistics SangerLCMFilteringSingularity_latest.sif -a filtered.annovar -b ${bam} -t $task.cpus > annotated.annovar
  singularity run --bind /nfs,/lustre --app additionalBAMStatistics SangerLCMFilteringSingularity_latest.sif -a annotated.annovar -b ${bam} -t $task.cpus -r $params.fasta \
    -s /lustre/scratch126/casm/team154pc/at31/chemo_trees/data/reference/snp_database > fully_annotated.annovar
  singularity run --bind /nfs,/lustre --app filtering SangerLCMFilteringSingularity_latest.sif -a fully_annotated.annovar -v ${caveman_vcf} -o out/ -p ${meta.donor_id}
  """
}

// Main workflow
workflow {

  // Print help message, supply typical command line usage for the pipeline
  if (params.help) {
    log.info paramsHelp("nextflow run nf-chemo-trees --sample_sheet sample_sheet.csv")
    exit 0
  }

  // validate the input parameters
  //validateParameters()

  // print summary of supplied parameters
  log.info paramsSummaryLog(workflow)

  // create a new channel of metadata from a sample sheet
  //ch_input = Channel.fromSamplesheet("sample_sheet")
  //ch_input.view()

  // get metadata + bam paths
  ch_input = \
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
          file(row.bam + ".met.gz", checkIfExists: true)
        ]
    }

  view(params.fasta)
  
  // run hairpin
  HAIRPIN(ch_input, params.fasta)

}