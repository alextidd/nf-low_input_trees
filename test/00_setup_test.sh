#!/bin/bash

# modules
module load samtools-1.19/python-3.12.0 

# files
ss=/lustre/scratch126/casm/team154pc/at31/chemo_trees/out/metadata/sample_sheet.csv

# TODO: add *.bas and *.met.gz files - how to generate?
# test dataset of chr22 subsamples from 3 samples from 1 donor
echo "donor_id,sample_id,bam,pindel_vcf,caveman_vcf" > test/sample_sheet.csv
while read -r donor_id sample_id bam pindel_vcf caveman_vcf ; do

  echo -e "donor: $donor_id\nsample: $sample_id"

  echo "subsampling BAM..."
  samtools view -b $bam chr22 > test/${sample_id}_chr22.bam
  samtools index $bam

  echo "subsampling CaVEMAN VCF..."
  zgrep -w '^#\|^#CHROM\|^chr22' $caveman_vcf | gzip \
  > test/${sample_id}_chr22.caveman.vcf.gz

  echo "subsampling Pindel VCF..."
  zgrep -w '^#\|^#CHROM\|^chr22' $pindel_vcf | gzip \
  > test/${sample_id}_chr22.pindel.vcf.gz

  echo "adding to the sample sheet..."
  echo "$donor_id,$sample_id,test/${sample_id}_chr22.bam,test/${sample_id}_chr22.pindel.vcf.gz,test/${sample_id}_chr22.caveman.vcf.gz" \
  >> test/sample_sheet.csv

  echo

done < <(cat $ss | head -4 | sed 1d | tr "," "\t")
