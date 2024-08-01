process sequoia_run {
  tag "${meta.donor_id}"
  label "week100gb"
  publishDir "${params.outdir}/${meta.donor_id}/", mode: 'copy'
  errorStrategy = 'retry'

  input:
  tuple val(meta),
        path(cgpVAF_out),
        path(fasta), path(fai)
  
  output:
    path("sequoia/*")

  script:
  """
  # plot trees of SNVs and indels separately
  build_phylogeny.R \
    --donor_id ${meta.donor_id} \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/ \
    --ncores ${task.cpus} \
    --beta_binom_shared ${params.sequoia_beta_binom_shared} \
    --normal_flt ${params.sequoia_normal_flt} \
    --snv_rho ${params.sequoia_snv_rho} \
    --indel_rho ${params.sequoia_indel_rho} \
    --min_cov ${params.sequoia_min_cov} \
    --max_cov ${params.sequoia_max_cov} \
    --only_snvs ${params.sequoia_only_snvs} \
    --keep_ancestral ${params.sequoia_keep_ancestral} \
    --exclude_samples ${params.sequoia_exclude_samples} \
    --cnv_samples ${params.sequoia_cnv_samples} \
    --vaf_absent ${params.sequoia_vaf_absent} \
    --vaf_present ${params.sequoia_vaf_present} \
    --mixmodel ${params.sequoia_mixmodel} \
    --min_clonal_mut ${params.sequoia_min_clonal_mut} \
    --tree_mut_pval ${params.sequoia_tree_mut_pval} \
    --genotype_conv_prob ${params.sequoia_genotype_conv_prob} \
    --min_pval_for_true_somatic ${params.sequoia_min_pval_for_true_somatic} \
    --min_variant_reads_shared ${params.sequoia_min_variant_reads_shared} \
    --min_vaf_shared ${params.sequoia_min_vaf_shared} \
    --create_multi_tree ${params.sequoia_create_multi_tree} \
    --germline_cutoff ${params.sequoia_germline_cutoff} \
    --plot_spectra ${params.sequoia_plot_spectra} \
    --max_muts_plot ${params.sequoia_max_muts_plot} \
    --lowVAF_filter ${params.sequoia_lowVAF_filter} \
    --lowVAF_filter_positive_samples ${params.sequoia_lowVAF_filter_positive_samples} \
    --VAF_treshold_mixmodel ${params.sequoia_VAF_treshold_mixmodel}
  
  # plot tree integrating SNVs and indels
  build_phylogeny.R \
    --donor_id ${meta.donor_id} \
    --cgpvaf_output ${cgpVAF_out.join(',')} \
    --genomeFile ${fasta} \
    --output_dir sequoia/ \
    --split_trees F \
    --ncores ${task.cpus} \
    --beta_binom_shared ${params.sequoia_beta_binom_shared} \
    --normal_flt ${params.sequoia_normal_flt} \
    --snv_rho ${params.sequoia_snv_rho} \
    --indel_rho ${params.sequoia_indel_rho} \
    --min_cov ${params.sequoia_min_cov} \
    --max_cov ${params.sequoia_max_cov} \
    --only_snvs ${params.sequoia_only_snvs} \
    --keep_ancestral ${params.sequoia_keep_ancestral} \
    --exclude_samples ${params.sequoia_exclude_samples} \
    --cnv_samples ${params.sequoia_cnv_samples} \
    --vaf_absent ${params.sequoia_vaf_absent} \
    --vaf_present ${params.sequoia_vaf_present} \
    --mixmodel ${params.sequoia_mixmodel} \
    --min_clonal_mut ${params.sequoia_min_clonal_mut} \
    --tree_mut_pval ${params.sequoia_tree_mut_pval} \
    --genotype_conv_prob ${params.sequoia_genotype_conv_prob} \
    --min_pval_for_true_somatic ${params.sequoia_min_pval_for_true_somatic} \
    --min_variant_reads_shared ${params.sequoia_min_variant_reads_shared} \
    --min_vaf_shared ${params.sequoia_min_vaf_shared} \
    --create_multi_tree ${params.sequoia_create_multi_tree} \
    --germline_cutoff ${params.sequoia_germline_cutoff} \
    --plot_spectra ${params.sequoia_plot_spectra} \
    --max_muts_plot ${params.sequoia_max_muts_plot} \
    --lowVAF_filter ${params.sequoia_lowVAF_filter} \
    --lowVAF_filter_positive_samples ${params.sequoia_lowVAF_filter_positive_samples} \
    --VAF_treshold_mixmodel ${params.sequoia_VAF_treshold_mixmodel}
  """
}

workflow sequoia {
  take:
  ch_input
  ch_fasta

  main:
  ch_input 
  | combine(ch_fasta)
  | sequoia_run
}