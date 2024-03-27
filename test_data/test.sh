#!/bin/bash
module add nextflow

# Run the workflow on the test data, and write the output to output/
nextflow-23.10.1-all \
    run \
    ../main.nf \
    --sample_sheet sample_sheet.csv \
    --genome_fasta genome_fasta/NC_001422.1.fasta \
    --outdir output \
    -with-report \
    -resume 