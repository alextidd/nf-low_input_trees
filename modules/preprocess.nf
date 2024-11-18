include { samplesheetToList } from 'plugin/nf-schema'

workflow preprocess {
  main:
  
  // TODO: change these reference files from channels to files
  // get fasta + fai
  fasta = \
    [file(params.fasta, checkIfExists: true),
     file(params.fasta + ".fai", checkIfExists: true)]

  // get bed + tbi
  cgpvaf_high_depth_bed = \
    [file(params.cgpvaf_high_depth_bed, checkIfExists: true),
     file(params.cgpvaf_high_depth_bed + ".tbi", checkIfExists: true)]

  // get cgpVAF normal bam + bai
  cgpvaf_normal_bam = \
    [file(params.cgpvaf_normal_bam),
     file(params.cgpvaf_normal_bam + ".bas", checkIfExists: true),
     file(params.cgpvaf_normal_bam + ".bai", checkIfExists: true),
     file(params.cgpvaf_normal_bam + ".met.gz", checkIfExists: true)]

  // get metadata + bam paths
  // TODO: check if we actually need bas and met files
  Channel.fromList(samplesheetToList(params.sample_sheet, "./assets/schema_sample_sheet.json"))
  | map { meta, bam, caveman_vcf, pindel_vcf ->
          [meta, 
          file(caveman_vcf),
          file(pindel_vcf),
          file(bam),
          // get bams' index, bas and met files
          file(bam + ".bai", checkIfExists: true),
          file(bam + ".bas", checkIfExists: true),
          file(bam + ".met.gz", checkIfExists: true)]}
  // get number of samples per donor
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta.subMap("donor_id"), meta.sample_id,
           caveman_vcf, pindel_vcf, bam, bai, bas, met] }
  | groupTuple
  | map { meta, sample_id, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta + [count: sample_id.size()], sample_id,
           caveman_vcf, pindel_vcf, bam, bai, bas, met] }
  | transpose
  | map { meta, sample_id, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta + [sample_id: sample_id],
           caveman_vcf, pindel_vcf, bam, bai, bas, met] }
  | set { ch_samples }

  // get caveman vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta + [vcf_type: "caveman"], caveman_vcf,
           bam, bai, bas, met]
  }
  | set { ch_caveman }

  // get pindel vcfs
  ch_samples 
  | map { meta, caveman_vcf, pindel_vcf, bam, bai, bas, met ->
          [meta + [vcf_type: "pindel"], pindel_vcf,
           bam, bai, bas, met]
  }
  | set { ch_pindel }

  // concat channels
  ch_caveman.concat(ch_pindel)
  | set { ch_input }

  emit:
  ch_input = ch_input
  fasta = fasta
  cgpvaf_high_depth_bed = cgpvaf_high_depth_bed
  cgpvaf_normal_bam = cgpvaf_normal_bam
}