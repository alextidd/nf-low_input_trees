process cgpVAF_run {
  tag "${meta.donor_id}:${vcf_type}:${chromosome}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"

  input: 
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals)
  tuple path(fasta), 
        path(fai)
  tuple path(high_depth_bed), 
        path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)
  each chromosome

  output:
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals)
  tuple path(fasta), 
        path(fai)
  tuple path(high_depth_bed), 
        path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)
  path "tmpvaf/*"

  script:
  def variant_type = (vcf_type == "caveman") ? "snp" : (vcf_type == "pindel") ? "indel" : ""
  """
  module load cgpVAFcommand/2.5.0

  # run cgpVAF
  cgpVaf.pl \
    --inputDir ./ \
    --outDir ./ \
    --variant_type ${variant_type} \
    --genome ${fasta} \
    --high_depth_bed ${params.high_depth_bed} \
    --bed_only 1 \
    --bedIntervals ${bed_intervals} \
    --base_quality 25 \
    --map_quality 30 \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.join(' ')} \
    -chr ${chromosome}
  
  # move to tmp
  mkdir tmpvaf/
  cp tmpvaf*/* tmpvaf/
  """
}

process cgpVAF_concat {
  tag "${meta.donor_id}:${vcf_type}:${chromosome}"
  label "normal"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"

  input: 
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals)
  tuple path(fasta), 
        path(fai)
  tuple path(high_depth_bed), 
        path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  script:  
  def variant_type = (vcf_type == "caveman") ? "snp" : (vcf_type == "pindel") ? "indel" : ""
  def ext = (vcf_type == "caveman") ? ".caveman_c.annot.vcf.gz" : (vcf_type == "pindel") ? ".pindel.annot.vcf.gz" : ""
  """
  cgpVaf.pl \
    --inputDir tmp/ \
    --outDir ./ \
    --variant_type ${variant_type} \
    --genome ${fasta} \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    -ct 1

  # rename output and remove vcf header
  cgpVAF_out=(`ls *_vaf.tsv`)
  n_files=\${#cgpVAF_out[@]}
  if (( n_files > 1 )) ; then
    echo "Error: More than one file matches pattern *_vaf.tsv"
    exit 1
  elif (( n_files == 1 )) ; then
    grep -v '^##' \${cgpVAF_out[0]} \
    > ${meta.donor_id}_${vcf_type}_vaf.tsv
  else
    echo "Error: No files match the pattern *_vaf.tsv"
    exit 1
  fi 
  """
}

workflow cgpVAF {
  take:
  ch_input
  ch_fasta
  ch_high_depth_bed
  ch_cgpVAF_normal_bam

  main:
  // create value channel of the chromosomes
  Channel.of(1..22).concat(Channel.of('X', 'Y')) | set { chromosomes }
  if (params.genome_build == 'hg38') {
   chromosomes = chromosomes | map { 'chr' + it } 
  }

  //run
  cgpVAF_run(
    ch_input,
    ch_fasta,
    ch_high_depth_bed,
    ch_cgpVAF_normal_bam,
    chromosomes)
  
  cgpVAF_run.out.view()

  emit:
  cgpVAF_run.out
}