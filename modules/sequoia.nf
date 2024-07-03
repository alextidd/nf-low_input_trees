process sequoia_run {
  tag "${meta.donor_id}"
  label "week100gb"
  publishDir "${params.outdir}/${meta.donor_id}/", mode: 'copy'
  errorStrategy = 'retry'

  input:
  tuple val(meta),
        path(cgpVAF_out)
  tuple path(fasta), path(fai)
  
  output:
    path("sequoia/*")

  script:
  """
  # plot trees of SNVs and indels separately
  Rscript build_phylogeny.R \
    --donor_id ${meta.donor_id} \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/ \
    --ncores $task.cpus
  
  # plot tree integrating SNVs and indels
  Rscript build_phylogeny.R \
    --donor_id ${meta.donor_id} \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/ \
    --split_trees F \
    --ncores $task.cpus
  # --plot_spectra T \
  """
}

workflow sequoia {
  take:
  ch_input
  ch_fasta

  main:
  sequoia_run(ch_input,
              ch_fasta)
}