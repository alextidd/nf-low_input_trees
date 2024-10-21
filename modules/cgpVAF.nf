process cgpVAF_run {
  tag "${meta.donor_id}:${meta.vcf_type}:${chr}"
  label "week50gb"
  maxRetries 3

  input: 
  tuple val(meta),
        val(sample_ids),
        path(vcfs), 
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals),
        val(chr)
  tuple path(fasta), path(fai)
  tuple path(high_depth_bed), path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam),
        path(cgpVAF_normal_bas),
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),
        path("tmpvaf_*/${chr}_progress.out"),
        path("tmpvaf_*/tmp_${chr}.tsv"),
        path("tmpvaf_*/tmp_${chr}.vcf")

  script:
  def variant_type = (meta.vcf_type == "caveman") ? "snp" : (meta.vcf_type == "pindel") ? "indel" : ""
  """
  module load cgpVAFcommand/2.5.0

  # gzip the vcfs
  for file in *.vcf ; do
    gzip -f \${file}
  done

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
    --vcf ${vcfs.collect { it + '.gz' }.join(' ')} \
    -chr ${chr}
  """
}

process cgpVAF_concat {
  tag "${meta.donor_id}:${meta.vcf_type}"
  label "week100gb"
  errorStrategy = "retry"
  publishDir "${params.outdir}/${meta.donor_id}/${meta.vcf_type}/", 
    mode: "copy"

  input: 
  tuple val(meta),
        val(sample_ids),
        path(vcfs),
        path(bams), path(bais), path(bass), path(mets),
        path(bed_intervals),
        path(tmpvaf_progress),
        path(tmpvaf_tsv),
        path(tmpvaf_vcf)
  tuple path(fasta), path(fai)
  tuple path(high_depth_bed), path(high_depth_tbi)
  tuple path(cgpVAF_normal_bam),
        path(cgpVAF_normal_bas),
        path(cgpVAF_normal_bai),
        path(cgpVAF_normal_met)

  output:
  tuple val(meta),   
        path("${meta.donor_id}_${meta.vcf_type}_vaf.tsv")
  
  script:
  def variant_type = (meta.vcf_type == "caveman") ? "snp" : (meta.vcf_type == "pindel") ? "indel" : ""
  """
  module load cgpVAFcommand/2.5.0

  # stage in tmpvaf dir 
  # (cgpVAF uses the first sample id as the output directory suffix)
  mkdir -p tmpvaf_${sample_ids[0]}
  mv tmp*{.vcf,.tsv} *_progress.out tmpvaf_${sample_ids[0]}

  # gzip the vcfs
  for file in *.vcf ; do
    gzip -f \${file}
  done

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
    --vcf ${vcfs.collect { it + '.gz' }.join(' ')} \
    -ct 1
  
  # rename output and remove vcf header
  cgpVAF_out=(`ls *_vaf.tsv`)
  n_files=\${#cgpVAF_out[@]}
  if (( n_files > 1 )) ; then
    echo "Error: More than one file matches pattern *_vaf.tsv"
    exit 1
  elif (( n_files == 1 )) ; then
    grep -v '^##' \${cgpVAF_out[0]} \
    > ${meta.donor_id}_${meta.vcf_type}_vaf.tsv
  else
    echo "Error: No files match the pattern *_vaf.tsv"
    exit 1
  fi 
  """
}

workflow cgpVAF {
  take:
  ch_input
  fasta
  high_depth_bed
  cgpVAF_normal_bam

  main:
  // create value channel of the chromosomes
  Channel.of(1..22).concat(Channel.of('X', 'Y')) | set { chromosomes }
  if (params.genome_build == 'hg38') {
   chromosomes = chromosomes | map { 'chr' + it }
  }

  // run cgpVAF by chromosome
  // TODO: figure out how to run cgpVAF in batches within patients
  cgpVAF_run(ch_input.combine(chromosomes),
             fasta, high_depth_bed, cgpVAF_normal_bam)

  // combine chromosomal channels
  cgpVAF_run.out
  | map { meta, chr_progress, chr_tsv, chr_vcf ->
          [groupKey(meta, 24), chr_progress, chr_tsv, chr_vcf] } 
  | groupTuple
  | set { ch_cgpVAF_chrs }

  // concat cgpVAF outputs, group channels by donor
  cgpVAF_concat(ch_input.combine(ch_cgpVAF_chrs, by: 0),
                fasta, high_depth_bed, cgpVAF_normal_bam)
  | map { meta, cgpVAF_out -> 
          [meta.subMap("donor_id"), cgpVAF_out]}
  | groupTuple
  | view
  | set { ch_cgpVAF_out }

  emit:
  ch_cgpVAF_out
}