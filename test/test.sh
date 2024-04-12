#!/bin/bash
wd=/lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/
cd $wd

# run the workflow on the test data, and write the output to out/
nextflow run ../main.nf \
    --sample_sheet data/sample_sheet.csv \
    --outdir out/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config
