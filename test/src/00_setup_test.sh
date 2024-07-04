#!/bin/bash

# use first 3 samples from 1st donor from sample sheet for test run
cat /lustre/scratch126/casm/team154pc/at31/chemo_trees/out/metadata/sample_sheet.csv |
head -4 \
> test/data/sample_sheet.csv