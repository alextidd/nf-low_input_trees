process sequoia_run {
  tag "${meta.donor_id}"
  label "normal"
  publishDir "${params.outdir}/phylogeny", mode: 'copy'
  conda '/nfs/users/nfs_a/at31/miniforge3/envs/sequoia'

  input:
  tuple val(meta),
        path(cgpVAF_out)
  tuple path(fasta), 
        path(fai)

  script:
  """
  build_phylogeny.R \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta}
  """
}

workflow sequoia {
  take:
  ch_input

  main:
  sequoia_run(
    ch_input)
}