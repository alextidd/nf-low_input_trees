#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/; ~/bin/jsub lsf -q week -n test_chemo_trees -m 2g -l log "bash test.sh" | bsub

# dirs
cd /lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/

module load singularity
# run the workflow on the hg38 test data
nextflow run ../main.nf \
    --sample_sheet data/GRCh38/chr22/sample_sheet.csv \
    --outdir out/GRCh38/chr22/ \
    -w work/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
    -resume

nextflow run ../main.nf \
    --sample_sheet data/GRCh38/chr10/sample_sheet.csv \
    --outdir out/GRCh38/chr10/ \
    -w work/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
    -resume
