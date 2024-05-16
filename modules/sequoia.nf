process sequoia_run {
  tag "${meta.donor_id}"
  label "normal4core"
  publishDir "${params.outdir}/${meta.donor_id}/", mode: 'copy'
  conda '/nfs/users/nfs_a/at31/miniforge3/envs/sequoia'

  input:
  tuple val(meta),
        path(cgpVAF_out)
  tuple path(fasta), path(fai)
  
  output:
    path("sequoia/*")

  script:
  """
  # plot trees of SNVs and indels separately
  build_phylogeny.R \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/
  
  # plot tree integrating SNVs and indels
  build_phylogeny.R \
    --donor_id ${meta.donor_id} \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/ \
    --split_trees F \
    --ncores $task.cpus
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