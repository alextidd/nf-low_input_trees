#!/bin/bash
# cd /lustre/scratch126/casm/team154pc/at31/low_input_trees; bsub -q week -M2000 -R "span[hosts=1] select[mem>2000] rusage[mem=2000]" -J low_input_trees -o test/log/low_input_trees_%J.out -e test/log/low_input_trees_%J.err "bash test/src/01_run_test.sh"

module load singularity

# run the nextflow test
nextflow run . \
  -profile test,sanger_hg38 \
  -w test/work/ \
  -resume