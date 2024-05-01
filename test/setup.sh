#!/bin/bash
wd=/lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/
mkdir -p data/GRCh{37,38}/chr22/
cd $wd

# modules
module load samtools-1.19.2/python-3.11.6

# generate GRCh38 test dataset
outdir=$wd/data/GRCh38/chr22/
head -1 ../../data/sample_sheet.csv > $outdir/sample_sheet.csv
while IFS=, read -r donor_id sample_id project_number bam pindel_vcf caveman_vcf ; do

    # subset vcfs to chr22
    new_vcfs=()
    for vcf in $pindel_vcf $caveman_vcf; do
        base_vcf=$(basename $vcf)
        new_vcf="$outdir/${base_vcf}"
        zcat $vcf | grep '^#\|^chr22' | gzip \
        > $new_vcf
        new_vcfs+=($new_vcf)
    done

    # subset bam to chr22
    base_bam=$(basename $bam)
    new_bam="$outdir/${base_bam}"
    samtools view -bo $new_bam $bam chr22

    # create bai, bas and met.gz files
    samtools index $new_bam
    samtools idxstats $new_bam > $new_bam.bas
    samtools flagstat $new_bam | gzip > $new_bam.met.gz

    # add to sample sheet
    echo "$donor_id,$sample_id,$project_number,$new_bam,${new_vcfs[0]},${new_vcfs[1]}" \
    >> $outdir/sample_sheet.csv

done < <(cat ../../data/sample_sheet.csv | sed 1d | head -2)

# generate GRCh37 test dataset
