#!bin/bash
set -euo pipefail

INPUT_DIR="results/fragment_lengths"
OUTPUT_DIR="results/fragment_stats"
mkdir -p "$OUTPUT_DIR"


for file in "$INPUT_DIR"/*_length_freq.tsv; do
    sample=$(basename "${file%_length_freq.tsv}")
    echo "Computing stats for $sample ..."

    awk '
    {
        n++;
        sum += $1;
        sumsq += $1 * $1;
        if (min=="" || $1 < min) min=$1;
        if (max=="" || $1 > max) max=$1;
        lengths[n]=$1;
    }
    END {
        # compute mean, stdev
        mean = sum / n;
        sd = sqrt(sumsq / n - mean * mean);

        # compute median
        asort(lengths);
        if (n % 2) median = lengths[(n+1)/2];
        else median = (lengths[n/2] + lengths[n/2+1]) / 2;

        printf "sample\tN\tmean\tmedian\tsd\tmin\tmax\n";
        printf "%s\t%d\t%.2f\t%.2f\t%.2f\t%d\t%d\n", sample, n, mean, median, sd, min, max;
    }' sample="$sample" "$file" > "$OUTPUT_DIR/${sample}_stats.tsv"
done
