#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow \
    run \
    ../main.nf \
    --sample_sheet sample_sheet.csv \
    --outdir out \
    --genome GRCh38 \
    -resume 