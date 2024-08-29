#!/bin/bash

# use 4 samples from 2 donors from sample sheet for test run (must test both
# inter-sample handling and inter-donor handling of channel behaviour)
# TODO: try to subsample (once met/bas requirement sorted out) for quicker tests
ss=/lustre/scratch126/casm/team154pc/at31/chemo_trees/out/metadata/sample_sheet.csv
( head -1 $ss ;
  grep 'PD63266ar\|PD63266ac\|PD63268bf\|PD63268ac' $ss ;
) | cat > test/data/sample_sheet.csv

# use 10 samples from 1 donor from sample sheet for test run
# (must test batching of samples in the cgpVAF step)
ss=/lustre/scratch126/casm/team154pc/at31/chemo_trees/out/metadata/sample_sheet.csv
( head -1 $ss ;
  grep 'PD63266' $ss | head -10 ;
) | cat > test/data/sample_sheet_10samples.csv