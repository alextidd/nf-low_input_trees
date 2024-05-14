#!/bin/bash
wd=/lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/
cd $wd

# modules
module load samtools-1.19.2/python-3.11.6

# generate GRCh38 test datasets (chr22 is small, chr10 is medium)
for chr in chr22 chrX ; do
    echo $chr
    outdir=$wd/data/GRCh38/$chr/
    mkdir -p $outdir

    head -1 ../../data/sample_sheet.csv > $outdir/sample_sheet.csv
    while IFS=, read -r donor_id sample_id project_id experiment_id bam pindel_vcf caveman_vcf ; do

        # subset vcfs to chr
        new_vcfs=()
        for vcf in $pindel_vcf $caveman_vcf; do
            base_vcf=$(basename $vcf)
            new_vcf="$outdir/${base_vcf}"
            zcat $vcf | grep "^#\|^$chr" | gzip \
            > $new_vcf
            new_vcfs+=($new_vcf)
        done

        # subset bam to chr
        base_bam=$(basename $bam)
        new_bam="$outdir/${base_bam}"
        samtools view -bo $new_bam $bam $chr

        # create bai, bas and met.gz files
        samtools index $new_bam
        samtools idxstats $new_bam > $new_bam.bas
        samtools flagstat $new_bam | gzip > $new_bam.met.gz

        # add to sample sheet
        echo "$donor_id,$sample_id,$project_id,$experiment_id,$new_bam,${new_vcfs[0]},${new_vcfs[1]}" \
        >> $outdir/sample_sheet.csv

    done < <(cat ../../data/sample_sheet.csv | sed 1d | head -2)
done

# test dataset for two full samples
mkdir -p data/GRCh38/full/
cat ../../data/sample_sheet.csv |
head -5 \
> data/GRCh38/full/sample_sheet.csv 