process sequoia_run {
  tag "${meta.donor_id}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/", mode: 'copy'
  conda '/nfs/users/nfs_a/at31/miniforge3/envs/sequoia'

  input:
  tuple val(meta),
        path(cgpVAF_out)
  tuple path(fasta), 
        path(fai)
  
  output:
    path("sequoia/*")

  script:
  """
  build_phylogeny.R \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/
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