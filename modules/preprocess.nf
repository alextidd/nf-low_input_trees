workflow preprocess {
  main:
  // get fasta + fai
  Channel.fromPath(params.fasta, checkIfExists: true)
  | map { fasta -> [fasta, file(fasta + ".fai", checkIfExists: true)] } 
  | set { ch_fasta }

  // get bed + tbi
  Channel.fromPath(params.high_depth_bed, checkIfExists: true)
  | map { bed -> [bed, file(bed + ".tbi", checkIfExists: true)] }
  | set { ch_high_depth_bed }

  // get cgpVAF normal bam + bai
  Channel.of(
    [file(params.cgpVAF_normal_bam, checkIfExists: true), 
    file(params.cgpVAF_normal_bas, checkIfExists: true),
    file(params.cgpVAF_normal_bai, checkIfExists: true),
    file(params.cgpVAF_normal_met, checkIfExists: true)])
  | set { ch_cgpVAF_normal_bam }

  // get metadata + bam paths
  Channel.fromPath(params.sample_sheet, checkIfExists: true)
  | splitCsv(header: true)
  | map { row ->
      meta = [sample_id: row.sample_id, donor_id: row.donor_id, 
              project_id: row.project_id, experiment_id: row.experiment_id]
      [meta, 
      file(row.caveman_vcf, checkIfExists: true),
      file(row.pindel_vcf, checkIfExists: true),
      file(row.bam, checkIfExists: true),
      file(row.bam + ".bai", checkIfExists: true),
      file(row.bam + ".bas", checkIfExists: true),
      file(row.bam + ".met.gz", checkIfExists: true)]
  }
  | set { ch_samples }

  // get caveman vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta, "caveman", caveman_vcf, bam, bai, bas, met]
  }
  | set { ch_caveman }

  // get pindel vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta, "pindel", pindel_vcf, bam, bai, bas, met]
  }
  | set { ch_pindel }

  // concatenate channels together
  ch_caveman.concat(ch_pindel) 
  | set { ch_input } 

  emit:
  ch_input = ch_input
  ch_fasta = ch_fasta
  ch_high_depth_bed = ch_high_depth_bed
  ch_cgpVAF_normal_bam = ch_cgpVAF_normal_bam
}