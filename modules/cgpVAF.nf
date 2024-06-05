process cgpVAF_run {
  tag "${meta.donor_id}:${vcf_type}:${chr}"
  label "week100gb"
  errorStrategy = 'retry'

  input: 
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals),
        path(fasta), path(fai),
        path(high_depth_bed), path(high_depth_tbi),
        path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met),
        val(chr)

  output:
  tuple val(meta),
        val(vcf_type),
        path("tmpvaf_*/${chr}_progress.out"),
        path("tmpvaf_*/tmp_${chr}.tsv"),
        path("tmpvaf_*/tmp_${chr}.vcf")

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
    --high_depth_bed ${high_depth_bed} \
    --bed_only 1 \
    --bedIntervals ${bed_intervals} \
    --base_quality 25 \
    --map_quality 30 \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.join(' ')} \
    -chr ${chr}
  """
}

process cgpVAF_concat {
  tag "${meta.donor_id}:${vcf_type}"
  label "week100gb"
  errorStrategy = "retry"
  publishDir "${params.outdir}/${meta.donor_id}/${vcf_type}/", 
    mode: "copy"

  input: 
  tuple val(meta),
        val(vcf_type), 
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals),
        path(tmpvaf_progress),
        path(tmpvaf_tsv),
        path(tmpvaf_vcf),
        path(fasta), path(fai),
        path(high_depth_bed), path(high_depth_tbi),
        path(cgpVAF_normal_bam), 
        path(cgpVAF_normal_bas), 
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),   
        path("${meta.donor_id}_${vcf_type}_vaf.tsv")
  
  script:
  def variant_type = (vcf_type == "caveman") ? "snp" : (vcf_type == "pindel") ? "indel" : ""
  """
  module load cgpVAFcommand/2.5.0
  
  # stage in tmpvaf dir 
  # (cgpVAF uses the first sample id as the output directory suffix)
  mkdir -p tmpvaf_${sample_ids[0]}
  mv tmp*{.vcf,.tsv} *_progress.out tmpvaf_${sample_ids[0]}

  # run cgpVAF concat
  cgpVaf.pl \
    --inputDir ./ \
    --outDir ./ \
    --variant_type ${variant_type} \
    --genome ${fasta} \
    --high_depth_bed ${high_depth_bed} \
    --bed_only 1 \
    --bedIntervals ${bed_intervals} \
    --base_quality 25 \
    --map_quality 30 \
    --tumour_name ${sample_ids.join(' ')} \
    --tumour_bam ${bams.join(' ')} \
    --normal_name "normal" \
    --normal_bam ${cgpVAF_normal_bam} \
    --vcf ${vcfs.join(' ')} \
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

  // run cgpVAF by chromosome
  ch_input 
  | combine(ch_fasta)
  | combine(ch_high_depth_bed)
  | combine(ch_cgpVAF_normal_bam)
  | combine(chromosomes) 
  | cgpVAF_run

  // combine chromosomal channels
  cgpVAF_run.out
  | groupTuple(by: [0,1])
  | set { ch_cgpVAF_chrs }

  // concat cgpVAF outputs
  ch_input
  | combine(ch_cgpVAF_chrs, by: [0,1])
  | combine(ch_fasta)
  | combine(ch_high_depth_bed)
  | combine(ch_cgpVAF_normal_bam)
  | cgpVAF_concat
  | groupTuple
  | set { ch_cgpVAF_out }

  emit:
  ch_cgpVAF_out
}