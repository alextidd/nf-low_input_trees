#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/low_input_sanger_filtering/test/; bsub -q week -M 2000 -R 'select[mem>=2000] span[hosts=1] rusage[mem=2000]' -J test_lisf -o log "bash test.sh"

# modules
module load singularity

# run the workflow on the hg38 test data, 3 full samples
nextflow run ../main.nf \
    --sample_sheet data/GRCh38/sample_sheet.csv \
    --project_type WES \
    --outdir out/GRCh38/ \
    -w ./work/ \
    -c ../conf/LSF.config \
    -resume \
    -N at31@sanger.ac.uk
