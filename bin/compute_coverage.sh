#!/bin/bash

bed=""
regions=""
OUTDIR=""
SAMPLE=""


# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --fragments) bed="$2"; shift ;;
        --chrom-sizes) chromSizes="$2"; shift ;;
        --regions) regions="$2"; shift ;;
        --sample) SAMPLE="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        --upstream-bp) upstream_bp="$2"; shift ;;
        --downstream-bp) downstream_bp="$2"; shift ;;
        --len-freq) len_freq="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

len_freq=results/${SAMPLE}/${SAMPLE}_length_freq.tsv
OUTDIR=results/${SAMPLE}
# echo bed is $bed 
# echo outdir is $OUTDIR
# echo sample is $SAMPLE
# echo regions is $regions

# Check required arguments
if [[  -z "$bed" || -z "$OUTDIR" || -z "$SAMPLE" || -z "$regions" ]]; then
    echo "Error: --fragments , --outdir , --regions and --sample are required."
    exit 1
fi

# sort regions and bed files
SORTED_REGIONS=$(echo "$regions" | sed 's/....$/.sorted.bed/')
fragments=$(echo "$bed" | sed 's/....$/.sorted.bed/')

echo sorting Regions and Fragments file  
LC_COLLATE=C sort -k1,1 -k2,2n $regions > $SORTED_REGIONS
LC_COLLATE=C sort -k1,1 -k2,2n $bed > $fragments
echo sorted regions is $SORTED_REGIONS

#get stats for the sample (mean, median, stdev)
echo running scripts/stats.R 
Rscript scripts/stats.R $len_freq results/$SAMPLE $SAMPLE

# confirm periodicity
lengths=$(echo "$bed" | sed 's/....$/_lengths.txt/')
echo running scripts/periodicity_confirmation.R
Rscript scripts/periodicity_confirmation.R $lengths $OUTDIR $SAMPLE


# echo "PWD        = $(pwd)"

#awk 'NR==FNR{next} {print "X", $0}'  test.txt
awk -v upstream_bp="$upstream_bp" \
-v downstream_bp="$downstream_bp" \
-v OUT_UP="data/regions_up_${upstream_bp}.bed" \
-v OUT_DOWN="data/regions_down_${downstream_bp}.bed" \
'BEGIN{OFS="\t"}
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

    # ----- upstream:  -----
    up_start = start - upstream_bp
    if (up_start < 0) up_start = 0
    up_end = start
    if (up_end > L) up_end = L
    if (up_end < up_start) up_end = up_start  # safety
    $1 = chr; $2 = up_start; $3 = up_end
    print > OUT_UP


    # ----- downstream: [end, end+2500) -----
    down_start = end
    if (down_start < 0) down_start = 0
    down_end = end + downstream_bp
    if (down_end > L) down_end = L
    if (down_end < down_start) down_end = down_start  # safety
    $1 = chr; $2 = down_start; $3 = down_end
    print > OUT_DOWN
}' $chromSizes  $SORTED_REGIONS


#bedtools map
OUT_UP="data/regions_up_${upstream_bp}.bed"
OUT_DOWN="data/regions_down_${downstream_bp}.bed"
echo overlap with regions FOR $SAMPLE ...
bedtools map \
  -a $SORTED_REGIONS \
  -b $fragments \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTDIR/${SAMPLE}_regions.bed

#bedtools map upper and lower
bedtools map \
  -a $OUT_UP \
  -b $fragments \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTDIR/${SAMPLE}_${upstream_bp}_up.bed

bedtools map \
  -a $OUT_DOWN \
  -b $fragments \
  -c 5,5,5,5 \
  -o count,mean,median,stdev \
  > $OUTDIR/${SAMPLE}_regions_${downstream_bp}.bed

Up="$OUTDIR/${SAMPLE}_${upstream_bp}_up.bed"
Down="$OUTDIR/${SAMPLE}_regions_${downstream_bp}.bed"
# normalize to uppper and lower

paste $OUTDIR/${SAMPLE}_regions.bed $Up $Down \
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
  }' > $OUTDIR/${SAMPLE}_regions_normalized_cov.bed

# delete bam files
rm  $OUTDIR/${SAMPLE}.deduped.bam
rm  $OUTDIR/${SAMPLE}.deduped.bam.bai

echo generating plots of normalized coverages on $SAMPLE
Rscript scripts/04_plot_cov_norm.R $SAMPLE $OUTDIR

