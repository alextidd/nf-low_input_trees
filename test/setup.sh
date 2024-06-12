#!/bin/bash

# dirs
bd=/lustre/scratch126/casm/team154pc/at31/
cd $bd/low_input_sanger_filtering/

# test dataset for 3 full samples
mkdir -p test/data/GRCh38/
cat $bd/chemo_trees/out/metadata/sample_sheet.csv |
head -4 \
> test/data/GRCh38/sample_sheet.csv 