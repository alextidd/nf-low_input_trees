#!/bin/bash
wd=/lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/
mkdir -p data/GRCh{37,38}
cd $wd

# generate GRCh38 test dataset
head -1 ../../data/sample_sheet.csv > data/GRCh38/sample_sheet.csv
while IFS=, read -r donor_id sample_id project_number bam pindel_vcf caveman_vcf ; do

    # subset vcfs to chr22
    new_vcfs=()
    for vcf in $pindel_vcf $caveman_vcf; do
        base_vcf=$(basename $vcf)
        new_vcf="$wd/data/GRCh38/${base_vcf/$sample_id/$sample_id.chr22}"
        zcat $vcf | grep -w '^#\|^chr22' | gzip \
        > $new_vcf
        new_vcfs+=($new_vcf)
    done

    # subset bam to chr22
    base_bam=$(basename $bam)
    new_bam="$wd/data/GRCh38/${base_bam/$sample_id/$sample_id.chr22}"
    samtools view -bo $new_bam $bam chr22

    # create bai, bas and met.gz files
    samtools index $new_bam
    samtools idxstats $new_bam > $new_bam.bas
    samtools flagstat $new_bam | gzip > $new_bam.met.gz

    # add to sample sheet
    echo "$donor_id,$sample_id,$project_number,$new_bam,${new_vcfs[0]},${new_vcfs[1]}" \
    >> data/GRCh38/sample_sheet.csv

done < <(cat ../../data/sample_sheet.csv | sed 1d | head -2)

# generate GRCh37 test dataset
