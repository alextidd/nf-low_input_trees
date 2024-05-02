process get_mutation_filtering_parameters {
  tag "${meta.donor_id}:${vcf_type}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${run}", 
    mode: "copy"

  input:
  tuple val(meta),
        val(sample_ids),
        path(merged_caveman),
        path(merged_pindel)
        tuple val(run), val(extra_arguments)

  output:
  tuple val(meta),
        val(vcf_type),
        val(sample_ids)
  
  script:
  """
  get_mutation_filtering_parameters.R \
    --runid ${meta.experiment_id} \
    --snvfile ${merged_caveman} \
    --indelfile ${merged_pindel} \
    --output ./ ${extra_arguments}
  """
}

workflow get_mutation_filtering_parameters {
    // get mutation filtering parameters
    cgpVAF.out 
    | branch { meta, vcf_type, merged_muts ->
                caveman: vcf_type == "caveman"
                    return tuple(meta, merged_muts)
                pindel: vcf_type == "pindel"
                    return tuple(meta, merged_muts) }
    | set { ch_merged_muts }

    // run twice:
    // 1. not excluding any samples
    // 2. excluding samples with peak VAF < 0.4
    runs = [["include_all", ""],
            ["remove_below_0.4", "--mixed_remove --vaf_cutoff 0.4"]]
    ch_merged_muts.caveman
    | join(ch_merged_muts.pindel)
    | combine(runs)
    | view
}

