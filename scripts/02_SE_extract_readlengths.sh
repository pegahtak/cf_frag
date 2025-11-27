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

echo "Extracting read lengths for $SAMPLE ..."

bam=$INPUT_DIR/${SAMPLE}.deduped.bam
samtools view -F 3840 "$bam" | \
        awk '{print length($10)}' > "${OUTPUT_DIR}/${SAMPLE}_lengths.txt"

echo "Generating bed file for $SAMPLE ..."
samtools view -F 3840 -q 30 "$bam" | \
        gawk '{
            chr=$3;
            start=$4;
	    len=length($10);
            end=$4 + length($10);
            strand = and($2,16) ? "-" : "+";
            print chr "\t" start "\t" end "\t" NR "\t" len "\t" strand;
        }' > "${OUTPUT_DIR}/${SAMPLE}.bed"

done < $SAMPLES_FILE
