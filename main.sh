#!/bin/bash

CONFIG=${1:-config/config.yaml}
source scripts/utils.sh

# read arguments from config file
REF=$(yq -r '.arguments.reference' <  "$CONFIG")
INPUT_DIR=$(yq -r '.arguments.input_dir' < "$CONFIG")
OUTPUT_DIR=$(yq -r '.arguments.output_dir' < "$CONFIG")
SAMPLES_FILE=$(yq -r '.arguments.samples_file' < "$CONFIG")
THREADS=$(yq -r '.arguments.threads' < "$CONFIG" )
DEDUP=$(yq -r '.arguments.deduplicate' < "$CONFIG" )
PICARD=$(yq -r '.arguments.picard' < "$CONFIG" )
PAIRED=$(yq -r '.arguments.paired_assay' < "$CONFIG" )
REGIONS=$(yq -r '.arguments.regions' < "$CONFIG" )
TARGETED=$(yq -r '.arguments.targeted' < "$CONFIG" )
DATADIR=$(yq -r '.arguments.datadir' < "$CONFIG" )


mkdir -p results

#while read SAMPLE; do
# 
# #Add or Replace read Groups
# java -jar $PICARD AddOrReplaceReadGroups I=$INPUT_DIR/${SAMPLE}.bam \
# O=$INPUT_DIR/${SAMPLE}.rg.bam \
# RGID=group1 \
# RGLB=lib1 \
# RGPL=ILLUMINA \
# RGPU=unit1 \
# RGSM=$SAMPLE
# 
# if [ "$DEDUP" == "true" ]; then
# 
# bash scripts/01_deduplicate.sh --PICARD $PICARD \
# --INPUT $INPUT_DIR/${SAMPLE}.bam \
# --OUTPUT $OUTPUT_DIR/${SAMPLE}.deduped.bam --REF $REF --SAMPLE $SAMPLE
# 
# fi
# 
# 
# done < $SAMPLES_FILE
while read -r SAMPLE; do
bash scripts/01_deduplicate.sh --PICARD $PICARD \
--INPUT $INPUT_DIR/${SAMPLE}.bam \
--OUTPUT $INPUT_DIR/${SAMPLE}.deduped.bam  --SAMPLE $SAMPLE
done < $SAMPLES_FILE
# echo extract fragments 
# 
# if [ $PAIRED == "true" ]; then
# bash scripts/02_PE_extract_fragments.sh --SAMPLES_FILE "$SAMPLES_FILE"
# else
# bash scripts/02_SE_extract_readlengths.sh --SAMPLES_FILE "$SAMPLES_FILE"
# fi 
# 
# echo calculating Total fragment stats ...
# Rscript scripts/03_fragment_length_stats.R "results/fragment_lengths"
# 
# if [ $TARGETED == "true" ]; then
# echo calculating fragment stats on regions of interest ...
# 
# bash scripts/04_regions_stats.sh --SAMPLES_FILE "$SAMPLES_FILE" --ASSAY_TYPE $PAIRED --REGIONS $REGIONS
# 
# fi 
# # echo FFT on fragment lengths
# # Rscript scripts/04_periodicity.R "results/fragment_lengths"
# 
# #echo calculate WPS
# #make WPS results folder
# #mkdir results/WPS
# #bash scripts/05_WPS.sh --SAMPLES_FILE "$SAMPLES_FILE" --REGIONS "$REGIONS"
# #BatchWPSfromCoordinateFile.sh  -i results/case1.deduped.bam -c $REGIONS
# # done< $SAMPLES_FILE
