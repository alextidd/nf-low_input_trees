#!/bin/bash

# Run the workflow on the test data, and write the output to output/
nextflow run ../main.nf \
    --sample_sheet sample_sheet.csv \
    --outdir out \
    --fasta /lustre/scratch125/casm/team268im/al28/bed_ref/GRCh38_full_analysis_set_plus_decoy_hla_genome.fa \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config