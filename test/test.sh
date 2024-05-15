#!/bi`n`/bash
# cd /lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/; bsub -q week -M 2000 -R 'select[mem>=2000] span[hosts=1] rusage[mem=2000]' -J test_chemo_trees -o log "bash test.sh"

# dirs
cd /lustre/scratch126/casm/team154pc/at31/chemo_trees/nf-chemo-trees/test/

# modules
module load singularity

# # run the workflow on the hg38 test data, chr22
# nextflow run ../main.nf \
#    --sample_sheet data/GRCh38/chr22/sample_sheet.csv \
#    --outdir out/GRCh38/chr22/ \
#    -w work/ \
#    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
#    -resume

# run the workflow on the hg38 test data, two full samples
nextflow run ../main.nf \
    --sample_sheet data/GRCh38/full/sample_sheet.csv \
    --outdir out/GRCh38/full/ \
    -w work/ \
    -c /lustre/scratch125/casm/team268im/at31/RA_som_mut/scomatic/config/LSF.config \
    -resume
