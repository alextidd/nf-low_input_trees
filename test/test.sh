#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/; mamba activate trees ; ~/bin/jsub lsf -q week -n test_chemo_trees -m 2g -l log "bash test.sh" | bsub

# dirs
wd=/lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/
cd $wd

# load singularity
module load singularity

# run the workflow on the hg38 test data
nextflow run ../main.nf \
    --sample_sheet data/GRCh38/sample_sheet.csv \
    --outdir out/GRCh38/ \
    -w work/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
    -resume

# run the workflow on hg37 test data
