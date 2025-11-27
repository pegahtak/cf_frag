#!/bin/bash 
set -euo pipefail



SAMPLES_FILE=""
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --SAMPLES_FILE) SAMPLES_FILE="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$SAMPLES_FILE" ]]; then
    echo "Error: --SAMPLES_FILE is required."
    exit 1
fi


INPUT_DIR="results/"
OUTPUT_DIR="results/fragment_lengths"
mkdir -p $OUTPUT_DIR

while IFS= read -r SAMPLE; do

echo "Extracting fragment lengths for $SAMPLE ..."

bam=$INPUT_DIR/${SAMPLE}.deduped.bam
samtools view -f 3 -F 3852 -q 30 $bam | \
        awk '{if($9>0) print $9}' > ${OUTPUT_DIR}/${SAMPLE}_lengths.txt

echo "Generating bed file for $SAMPLE ..."
samtools view -f 3 -F 3852 -q 30 $bam | \
    gawk '{
     if ($9>0){ 
            chr=$3;
	    start=$4;
	    end=$4+$9;
            if (and($2,16)) { strand = "-" } else { strand = "+" } 
	    print chr "\t" start "\t" end "\t" NR "\t" $9 "\t" ".";
       }
        
    }' > ${OUTPUT_DIR}/${SAMPLE}.bed

echo getting reads with lengths below 40bp ...

samtools view -h -f 3 -F 3852 "$bam" | \
awk 'BEGIN {OFS="\t"} 
     /^@/ {print; next} 
     ($9 > 0 && $9 < 40)' >  ${OUTPUT_DIR}/${SAMPLE}_short.sam

echo "Extracting fragment lengths for $SAMPLE ..."

Rscript scripts/fragment_lengths.R "${OUTPUT_DIR}/${SAMPLE}_lengths.txt" "${OUTPUT_DIR}/${SAMPLE}_length_freq.tsv"

done < $SAMPLES_FILE
