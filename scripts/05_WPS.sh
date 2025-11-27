#!/bin/bash 
set -euo pipefail



SAMPLES_FILE=""
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --SAMPLES_FILE) SAMPLES_FILE="$2"; shift ;;
        --REGIONS) REGIONS="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$SAMPLES_FILE" || -z $REGIONS ]]; then
    echo "Error: --SAMPLES_FILE and --REGIONS are required."
    exit 1
fi

INPUT_DIR="results"
OUTPUT_DIR="results/WPS"
mkdir -p $OUTPUT_DIR


while IFS= read -r SAMPLE; do

echo "Extracting fragment lengths for $SAMPLE ..."

samtools index ${INPUT_DIR}/$SAMPLE.deduped.bam
bash BatchWPSfromCoordinateFile.sh --average false  -i ${INPUT_DIR}/$SAMPLE.deduped.bam -c $REGIONS
echo Normalize WPS for ${SAMPLE}
bash CalcAverageWPS.sh --keep true \
-i  ${INPUT_DIR}/${SAMPLE}.deduped.bam_NF_TSS_coordinates_for_10k_genes+-5kb.txt_120WPS_1-500bp \
--num 1
done < $SAMPLES_FILE

for f in $(ls "${INPUT_DIR}" | grep 'WPS' | grep -v '^WPS$'); do
    mv "${INPUT_DIR}/${f}" "${OUTPUT_DIR}/"
done

#for f in $(ls ${OUTPUT_DIR}); do 
#awk '{sub(/^chrchr/, "chr", $1); print}' ${OUTPUT_DIR}/$f > ${OUTPUT_DIR}/$f.tmp
#mv ${OUTPUT_DIR}/$f.tmp ${OUTPUT_DIR}/$f
#done
