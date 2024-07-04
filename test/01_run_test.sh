#!/bin/bash

module load singularity

# run the nextflow test
nextflow run . -profile test --outdir out/test/