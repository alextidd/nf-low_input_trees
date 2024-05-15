#Load the required modules
module load cgpVAFcommand
module load mpboot
module load perl
module load vagrent
module load pcap-core
module load picard-tools
module load canPipe/live
module load canned-queries-client
module load bwa
module load samtools
module load hairpin

#Set primary variables
DONOR_ID=PD55782
LONG_DONOR_ID=PD55782_75M1
ALL_PROJECT_NUMBERS=2911 #If have more than one project number , can have as comma separated
EXP_ID=set1
DONOR_AGE=75 #Age is used at the end to create ultrametric trees scaled by age

#Secondary variables
PD_NUMBERS=$DONOR_ID #If have more than one PD number (e.g. for HSCT), can have as comma separate
MATS_AND_PARAMS_DIR=/lustre/scratch119/casm/team154pc/ld18/chemo/filtering_runs/mats_and_params
ROOT_DIR=/lustre/scratch119/casm/team154pc/ld18/chemo/${LONG_DONOR_ID}/${EXP_ID}
SNV_BEDFILE_NAME=chemo_${LONG_DONOR_ID}_${EXP_ID}_caveman.bed
INDEL_BEDFILE_NAME=chemo_${LONG_DONOR_ID}_${EXP_ID}_pindel.bed
MS_FILTERED_BEDFILE_NAME=chemo_${LONG_DONOR_ID}_${EXP_ID}_postMS_SNVs.bed

RUN_ID=$EXP_ID
RUN_ID_M=${EXP_ID}
RUN_ID_TB=${RUN_ID_M}
IFS=',' read -r -a PROJECT_ARRAY <<< "$ALL_PROJECT_NUMBERS"

#--------------------Submit LCM filtering jobs-----------------------------
# download QC excel file from Canapps & check if you're keeping all samples - note down if you're not!
# to get the complete.vcf files - use those files which end in  .caveman_c.annot.vcf.gz at /nfs/cancer_ref01/nst_links/live/${ALL_PROJECT_NUMBERS}/${PD_NUMBERS}*/
# go to local computer and prepare haripin commands files
# ~/Data/chemo/scripts/Sub_filter/generate_commands_for_subs_filtering_etc.R
# scp ~/Data/chemo/scripts/Sub_filter/commands/${DONOR_ID}_*flags.txt ld18@farm5-login:/lustre/scratch119/casm/team154pc/ld18/chemo/${LONG_DONOR_ID}/${EXP_ID}/
mkdir -p $ROOT_DIR/hairpin/output_files
cd $ROOT_DIR/
bsub -J"run[1-`wc -l<PD55782_75M1_v_flags.txt`]%50" -o cmd.%J.%I.out -e cmd.%J.%I.err -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" -M12000 -n 12 ./wrapper.sh PD55782_75M1_v_flags.txt PD55782_75M1_b_flags.txt PD55782_75M1_o_flags.txt

bsub -J"run[1-`wc -l<${LONG_DONOR_ID}_v_flags.txt`]%50" -o cmd.%J.%I.out -e cmd.%J.%I.err -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" -M12000 -n 12 ./wrapper.sh ${LONG_DONOR_ID}_v_flags.txt ${LONG_DONOR_ID}_b_flags.txt ${LONG_DONOR_ID}_o_flags.txt

# bsub -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" -M12000 -o hairpin/out.%J.log -e hairpin/err.%J.log -n 12 ./wrapper.sh ${LONG_DONOR_ID}_v_flags.txt ${LONG_DONOR_ID}_b_flags.txt ${LONG_DONOR_ID}_o_flags.txt

# bsub -R "select[mem>12000] rusage[mem=12000] span[hosts=1]" -M12000 -o hairpin/out.%J.log -e hairpin/err.%J.log -n 12 ../.././wrapper.sh PD55782_75M1_v_flags.txt PD55782_75M1_b_flags.txt PD55782_75M1_o_flags.txt
# /lustre/scratch119/casm/team154pc/ms56/my_programs/Submitting_Mathijs_filters_jobs.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -o $ROOT_DIR/MS_filters/output_files -q normal

# if files are generated in ${ROOT_DIR} move them to the hairpin directory

#--------------------SNV analysis-----------------------------
mkdir -p $ROOT_DIR/caveman_raw/caveman_pileup/output; cd $ROOT_DIR/caveman_raw
/lustre/scratch119/casm/team154pc/ms56/my_programs/import_new_samples_only.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t SNVs
cut -f 1,2,4,5 *.caveman_c.annot.vcf_pass_flags|sort|uniq>caveman_pileup/$SNV_BEDFILE_NAME
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Run createVafCmd - with input value '3' (this is the option for selecting the caveman files as input)
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "3" | 
    createVafCmd.pl \
        -pid $PROJECT_NUMBER  \
        -o output \
        -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
        -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz \
        -mq 30  -bo 1 \
        -b $SNV_BEDFILE_NAME
done
# echo "3"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/genome.fa -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz -mq 30  -bo 1  -b $SNV_BEDFILE_NAME

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/caveman_raw/caveman_pileup/output
# /lustre/scratch119/casm/team154pc/ms56/my_programs/create_split_config_ini.R -p $PROJECT_NUMBER
/lustre/scratch119/casm/team154pc/ld18/chemo/scripts/create_split_config_ini.R -p $PROJECT_NUMBER
cd $ROOT_DIR/caveman_raw/caveman_pileup

#Then re-run the createVafCmd.pl script with the new config file  - with input value '3' (this is the option for selecting the caveman files as input)
echo "3"| 
createVafCmd.pl \
    -pid $PROJECT_NUMBER  \
    -o output \
    -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini \
    -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz \
    -mq 30  -bo 1  \
    -b $SNV_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together and get more memory
sed -e 's/\%5/\%50/g;s/2000/4000/g;s/500/1000/g;' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#--------------------Indel analysis-----------------------------
mkdir -p $ROOT_DIR/pindel_raw/pindel_pileup/output; cd $ROOT_DIR/pindel_raw
/lustre/scratch119/casm/team154pc/ms56/my_programs/import_new_samples_only.R -p $ALL_PROJECT_NUMBERS -s $PD_NUMBERS -t indels
cut -f 1,2,4,5 *.pindel.annot.vcf_pass_flags|sort|uniq>pindel_pileup/$INDEL_BEDFILE_NAME
cd $ROOT_DIR/pindel_raw/pindel_pileup

#Run createVafCmd for each project id containing samples - with input value '1' (this is the option for selecting the pindel files as input)
for PROJECT_NUMBER in "${PROJECT_ARRAY[@]}"; do
    echo "1"|
    createVafCmd.pl \
    -pid $PROJECT_NUMBER \
    -o output \
    -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa \
    -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz \
    -mq 30  -bo 1  \
    -b $INDEL_BEDFILE_NAME
done
# echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/genome.fa -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME

#Then run the "create_split_config_ini.R" script in the output folder
PROJECT_NUMBER=${PROJECT_ARRAY[0]}
cd $ROOT_DIR/pindel_raw/pindel_pileup/output
/lustre/scratch119/casm/team154pc/ld18/chemo/scripts/create_split_config_ini.R -p $PROJECT_NUMBER
cd $ROOT_DIR/pindel_raw/pindel_pileup

#Then re-run the createVafCmd.pl script with the new config file - with input value '1' (this is the option for selecting the pindel files as input)
# echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/genome.fa -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh37d5/shared/ucscHiDepth_0.01_mrg1000_no_exon_coreChrs.bed.gz -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME
echo "1"|createVafCmd.pl -pid $PROJECT_NUMBER  -o output -i output/${PROJECT_NUMBER}_cgpVafConfig_split.ini -g /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/genome.fa -hdr /lustre/scratch119/casm/team78pipelines/reference/human/GRCh38_full_analysis_set_plus_decoy_hla/shared/HiDepth_mrg1000_no_exon_coreChrs_v3.bed.gz -mq 30  -bo 1  -b $INDEL_BEDFILE_NAME

#Update the run_bsub.sh command to allow more jobs in the array to run together & get more memory
sed -e 's/\%5/\%50/g;s/2000/4000/g' run_bsub.sh >run_bsub_updated.sh

#Run this if need to switch to the long queue (not normally necessary)
#sed -i -e 's/normal/long/g' run_bsub_updated.sh

bash run_bsub_updated.sh

#--------------------ONCE cgpVAF HAS COMPLETED-----------------------------
#-----------------------------SNV merge-----------------------------
# bsub -o $PWD/log.%J \
#     -e $PWD/err.%J \
#     -q normal \
#     -G team78-grp \
#     -R 'select[mem>=8000] span[hosts=1] rusage[mem=8000]' \
#     -M8000 \
#     -n1 \
#     -J "SNV_merge" \
#     /lustre/scratch119/casm/team154pc/ms56/my_programs/INDEL_merge.sh $EXP_ID $ROOT_DIR

cd $ROOT_DIR/caveman_raw/caveman_pileup/output/output/PDv38is_wgs/snp
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,24,26,39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $(sed -n '1p' files) > temp.1.cut   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do    
    if [ -s temp.$FILE ]
    then
        echo "temp file temp.$FILE already exists. Moving onto next file..."
    else
        echo "temp file temp.$FILE does not yet exist, will be created"
        cut -f 39,41,54,56,69,71,84,86,99,101,114,116,129,131,144,146,159,161,174,176 $FILE >  temp.$FILE.cut
    fi       
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    echo $FILE
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_SNVs_${EXP_ID}.tsv && rm output.*

mv $ROOT_DIR/caveman_raw/caveman_pileup/output/output/PDv38is_wgs/snp/merged_SNVs_${EXP_ID}.tsv $ROOT_DIR/

#-----------------------------INDEL merge-----------------------------
cd $ROOT_DIR/pindel_raw/pindel_pileup/output/output/PDv38is_wgs/indel
ls *_vaf.tsv > files

#for first file
cut -f 3,4,5,6,16,18,25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $(sed -n '1p' files) > temp.1   #Chr, Pos, Ref, Alt, + MTR, DEP for all samples (max 11 samples in one file)

#for subsequent files (exclude Chr, Pos, Ref, Alt and PDv37is)
for FILE in $(tail -n+2 files); do
    cut -f 25,27,34,36,43,45,52,54,61,63,70,72,79,81,88,90,97,99,106,108,115,117 $FILE >  temp.$FILE
done

#Remove empty rows where header was with awk
for FILE in $(ls temp.*); do
    awk 'NF' $FILE > output.$FILE
    rm $FILE
done

#Concatenate output files to one merged file & move to the root directory
paste output.* > merged_indels_${EXP_ID}.tsv && rm output.*
mv $ROOT_DIR/pindel_raw/pindel_pileup/output/output/PDv38is_wgs/indel/merged_indels_${EXP_ID}.tsv $ROOT_DIR/

#--------------------AFTER ALL CGPVAF MATRICES ARE GENERATED-----------------------------
cd $ROOT_DIR
mkdir -p log_files
mkdir -p err_files

#Run the "Mutation_filtering_get_params" script. Run twice: (1) not excluding samples, (2) excluding colonies with peak VAF <0.4
bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q basement -G team154-vwork -R 'select[mem>=32000] span[hosts=1] rusage[mem=32000]' \
    -M32000 -n6 -J GET_PARAMS \
    /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/Mutation_filtering_get_parameters.R \
    -r $RUN_ID \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR

bsub -o $ROOT_DIR/log_files/mats_and_params.log.%J -e $ROOT_DIR/err_files/mats_and_params.err.%J \
    -q basement -G team154-vwork -R 'select[mem>=32000] span[hosts=1] rusage[mem=32000]' \
    -M32000 -n6 -J GET_PARAMS \
    /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/Mutation_filtering_get_parameters.R  \
    -r $RUN_ID_M \
    -s $ROOT_DIR/merged_SNVs_${EXP_ID}.tsv \
    -i $ROOT_DIR/merged_indels_${EXP_ID}.tsv \
    -o $MATS_AND_PARAMS_DIR \
    -m \
    -v 0.4

#ONCE MS FILTERS JOBS HAVE FINISHED
# cd $ROOT_DIR/MS_filters/output_files
cd $ROOT_DIR/hairpin/
# perl /lustre/scratch119/casm/team154pc/ms56/my_programs/filter_for_bed_file.pl
perl /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/filter_for_bed_file.pl
# cut -f 1,2,4,5 *complete_final_retained_3.vcf_for_bed_file|sort|uniq>$ROOT_DIR/MS_filters/$MS_FILTERED_BEDFILE_NAME
cut -f 1,2,4,5 *complete.hairpin.vcf_for_bed_file|sort|uniq>$ROOT_DIR/hairpin/$MS_FILTERED_BEDFILE_NAME

#ONCE Mutation_filtering_get_paramaters.R SCRIPT HAS COMPLETED
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/Reducing_mutset_from_MSfilters.R -r $RUN_ID -b $ROOT_DIR/hairpin/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR
bsub -o $ROOT_DIR/log_files/MSReduce.log.%J -e $ROOT_DIR/err_files/MSReduce.err.%J -q yesterday -R 'select[mem>=16000] span[hosts=1] rusage[mem=16000]' -M16000 -n1 -J MSreduce /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/Reducing_mutset_from_MSfilters.R -r $RUN_ID_M -b $ROOT_DIR/hairpin/$MS_FILTERED_BEDFILE_NAME -d $MATS_AND_PARAMS_DIR

#Run the sensitivity analysis (if run on LCM pathway)

# bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
#    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
#    -M4000 -n1 -J sensitivity \
#    /lustre/scratch119/casm/team154pc/ms56/my_programs/Sensitivity_analysis_from_SNPs.R \
#    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID}_postMS \
#    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
#    -i $ROOT_DIR/pindel_raw \
#    -s $ROOT_DIR/MS_filters/output_files \
#    -x '_complete_final_retained_3.vcf_for_bed_file'

bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
    -q normal -R 'select[mem>=40000] span[hosts=1] rusage[mem=40000]' \
    -M40000 -n1 -J sensitivity \
    /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/Sensitivity_analysis_from_SNPs.R \
    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID}_postMS \
    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
    -i $ROOT_DIR/pindel_raw \
    -s $ROOT_DIR/hairpin/ \
    -x 'complete.hairpin.vcf_for_bed_file'


#Run the sensitivity analysis (if not on LCM pathway)
# bsub -o $ROOT_DIR/log_files/sensitivity.log.%J -e $ROOT_DIR/err_files/sensitivity.err.%J \
#    -q normal -R 'select[mem>=4000] span[hosts=1] rusage[mem=4000]' \
#    -M4000 -n1 -J sensitivity \
#    /lustre/scratch119/casm/team154pc/ms56/my_programs/Sensitivity_analysis_from_SNPs.R \
#    -m $MATS_AND_PARAMS_DIR/mats_and_params_${RUN_ID} \
#    -o $ROOT_DIR -n sensitivity_analysis_${EXP_ID} \
#    -i $ROOT_DIR/pindel_raw \
#    -s $ROOT_DIR/caveman_raw \
#    -x '.caveman_c.annot.vcf_pass_flags'

#-----------------------Run the tree-building script - this has lots of options

# -i The id for the tree-building run - will be included in the output file names.
# -m Path to mutation filtering output file
# -f Option to do p-value based filtering, or vaf based filtering (pval or vaf)
# -c Cut off to exclude samples with low coverage
# -o Folder for storing script output files. Folder & subfolders will be created if don't already exist.
# -s Path to the sensitivity dataframe
# -p Save trees as polytomous trees (not bifurcating trees)
# -t The age ('time') of the individual - for building the age-adjusted ultrametric tree
# -j Option to do initial tree-building with just the SNVs (i.e. don't' include indels)
# -a Option to keep an ancestral branch

# RUN_ID_TB=${EXP_ID}_m40_reduced
RUN_ID_TB=${EXP_ID}_reduced

#bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
#    -q basement -R 'select[mem>=24000] span[hosts=1] rusage[mem=24000]' \
#    -M24000 -n1 -J tree_build \
#    /lustre/scratch119/casm/team154pc/ms56/my_programs/filtering_from_table_mix_remove.R \
#    -i ${RUN_ID_TB}_a_j \
#    -m /lustre/scratch119/casm/team154pc/ms56/chemo_exposed/filtering_runs/mats_and_params/mats_and_params_${RUN_ID_TB} \
#    -f pval \
#    -d $DONOR_ID \
#    -c 4 \
#    -o /lustre/scratch119/casm/team154pc/ms56/chemo_exposed/filtering_runs \
#    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
#    -p \
#    -t $DONOR_AGE \
#    -j \
#    -a 

# must export vcf header into file
# for 75M1 this was
# ROOT_DIR:  cat ../old_sub/PD55782aa_complete.vcf | head -3416 > VCF_header_for_VaGrent.txt

# CAREFUL TO SET ThE MATS AND PARAMS FILE CORRECTLY

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q basement -R 'select[mem>=100000] span[hosts=1] rusage[mem=100000]' \
    -M100000 -n1 -J tree_build \
    /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/filtering_from_table_mix_remove.R \
    -i ${RUN_ID_TB}_a_j \
    -m /lustre/scratch119/casm/team154pc/ld18/chemo/filtering_runs/mats_and_params/mats_and_params_set1_postMS_reduced \
    -f pval \
    -d $DONOR_ID \
    -c 4 \
    -o $ROOT_DIR/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a

#bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J -q basement -R 'select[mem>=100000] span[hosts=1] rusage[mem=100000]' -M100000 -n1 -J tree_build /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/filtering_from_table_mix_remove.R -i ${RUN_ID_TB}_a_j -m /lustre/scratch119/casm/team154pc/ms56/chemo_exposed/filtering_runs/mats_and_params/mats_and_params_${RUN_ID_TB} -f pval -d $DONOR_ID -c 4 -o $ROOT_DIR/filtering_runs -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} -p -t $DONOR_AGE -j -a

bsub -o $ROOT_DIR/log_files/treebuild.log.%J -e $ROOT_DIR/err_files/treebuild.err.%J \
    -q basement -R 'select[mem>=150000] span[hosts=1] rusage[mem=150000]' \
    -M150000 -n1 -J tree_build \
    /lustre/scratch119/casm/team154pc/ld18/chemo/scripts/filtering_from_table_mix_remove.R \
    -i ${RUN_ID_TB} \
    -m /lustre/scratch119/casm/team154pc/ld18/chemo/filtering_runs/mats_and_params/mats_and_params_set1_postMS_reduced \
    -f vaf \
    -d $DONOR_ID \
    -c 4 \
    -o $ROOT_DIR/filtering_runs \
    -s $ROOT_DIR/sensitivity_analysis_${EXP_ID} \
    -p \
    -t $DONOR_AGE \
    -j \
    -a




