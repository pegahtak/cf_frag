#!/bin/bash
# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample-file) sample_file="$2"; shift ;;
         --regions) regions="$2"; shift ;;
         --upstream-bp) upstream_bp="$2"; shift ;;
         --downstream-bp) downstream_bp="$2"; shift ;;
         --chrom-sizes) chrom_sizes="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

mkdir -p batch_scripts

while IFS=$'\t' read -r sample group bam; do
#creat shebang
echo '#!/bin/bash' > batch_scripts/$sample.sh

echo python bin/cf_frags.py preprocess-bam --bam $bam --sample $sample \
--outdir results/$sample >> batch_scripts/$sample.sh
echo python bin/cf_frags.py coverage --fragments results/$sample/$sample.bed \
--regions $regions --sample $sample \
--upstream-bp $upstream_bp --downstream-bp $downstream_bp \
--outdir results/$sample --len-freq results/$sample/${sample}_len_freq.txt \
--chrom-sizes $chrom_sizes >> batch_scripts/$sample.sh

done< $sample_file

# run scripts in parallel
for s in batch_scripts/*.sh; do
  bash "$s" &
done
wait

python bin/cf_frags.py merge-plots --sample-file $sample_file --outdir results


