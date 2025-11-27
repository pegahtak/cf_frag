#!/bin/bash

set -euo pipefail


SAMPLES_FILE=""
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --SAMPLES_FILE) SAMPLES_FILE="$2"; shift ;;
        --ASSAY_TYPE) ASSAY_TYPE="$2"; shift;;
        --REGIONS) REGIONS="$2"; shift;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$SAMPLES_FILE" || -z $ASSAY_TYPE || -z $REGIONS ]]; then
    echo "Error: --SAMPLES_FILE, ASSAY_TYPE and REGIONS are required."
    exit 1
fi

INPUT_DIR="results/fragment_lengths"
OUTPUT_DIR="results/fragment_stats"


#sort regions
SORTED_REGIONS=$(echo "$REGIONS" | sed 's/....$/sorted.bed/')
echo $SORTED_REGIONS
echo sorting Regions file
LC_COLLATE=C sort -k1,1 -k2,2n $REGIONS > $SORTED_REGIONS
#LC_COLLATE=C sort -k1,1 -k2,2n $SORTED_REGIONS # check 

echo create 2.5kb window around regions
awk 'BEGIN{OFS="\t"}
NR==FNR {
    # Read chrom sizes
    chrlen[$1] = $2
    next
}
{
    chr   = $1
    start = $2
    end   = $3
    L     = chrlen[chr]

    # ----- upstream: [start-2500, start) -----
    up_start = start - 2500
    if (up_start < 0) up_start = 0
    up_end = start
    if (up_end > L) up_end = L
    if (up_end < up_start) up_end = up_start  # safety
    $1 = chr; $2 = up_start; $3 = up_end
    print > "data/regions_up_2500.bed"

    # ----- downstream: [end, end+2500) -----
    down_start = end
    if (down_start < 0) down_start = 0
    down_end = end + 2500
    if (down_end > L) down_end = L
    if (down_end < down_start) down_end = down_start  # safety
    $1 = chr; $2 = down_start; $3 = down_end
    print > "data/regions_down_2500.bed"
}' data/hg38.chrom.sizes  $SORTED_REGIONS

UPPER=data/regions_up_2500.bed
LOWER=data/regions_down_2500.bed

while IFS= read -r SAMPLE; do

# echo "Converting bam to bed for $SAMPLE ..."
#
# bam=$INPUT_DIR/${SAMPLE}.deduped.bam
bed=$INPUT_DIR/${SAMPLE}.bed
bedSort=$INPUT_DIR/${SAMPLE}.sorted.bed


#sort fragment file
echo sorting $SAMPLE bed file
LC_COLLATE=C sort -k1,1 -k2,2n $bed > $bedSort
# LC_COLLATE=C sort -k1,1 -k2,2n $bedSort

#bedtools map
echo overlap with regions FOR $SAMPLE ...
bedtools map \
  -a $SORTED_REGIONS \
  -b $bedSort \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTPUT_DIR/${SAMPLE}_regions.bed

#bedtools map upper and lower
bedtools map \
  -a $UPPER \
  -b $bedSort \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTPUT_DIR/${SAMPLE}_regions_up.bed

bedtools map \
  -a $LOWER \
  -b $bedSort \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTPUT_DIR/${SAMPLE}_regions_down.bed
  
# normalize to uppper and lower

paste $OUTPUT_DIR/${SAMPLE}_regions.bed $OUTPUT_DIR/${SAMPLE}_regions_up.bed $OUTPUT_DIR/${SAMPLE}_regions_down.bed \
  | awk 'BEGIN{OFS="\t"}{
      # Each file had 7 columns, so after paste:
      #  region:    $1..$7
      #  upstream:  $8..$14
      #  downstream:$15..$21

      cov_region = $7      # mean cov in region
      cov_up     = $17     # mean cov upstream
      cov_down   = $27     # mean cov downstream

      bg = (cov_up + cov_down) / 2.0

      if (bg > 0) {
          norm = cov_region / bg
      } else {
          norm = "NA"
      }

      # print region coordinates + raw and normalized values
      print $1, $2, $3, $4, $5, $6, cov_region, cov_up, cov_down, bg, norm
  }' > $OUTPUT_DIR/${SAMPLE}_regions_normalized_cov.bed
echo generating plots of normalized coverages on $SAMPLE
Rscript scripts/04_plot_cov_norm.R $SAMPLE

done < $SAMPLES_FILE
