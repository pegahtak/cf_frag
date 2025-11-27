#!/bin/bash

INPUT=""
OUTDIR=""
SAMPLE=""

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --bam) INPUT="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        --sample) SAMPLE="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[  -z "$INPUT" || -z "$OUTDIR" || -z "$SAMPLE" ]]; then
    echo "Error: --bam , --outdir and --sample are required."
    exit 1
fi
OUTDIR=results/${SAMPLE}
mkdir -p $OUTDIR
base=$OUTDIR/$SAMPLE

#AddOrReplaceReadGroups
RG_OUTPUT=$base.rg.bam
java -jar picard.jar AddOrReplaceReadGroups I=$INPUT \
O=$RG_OUTPUT \
RGID=group1 \
RGLB=lib1 \
RGPL=ILLUMINA \
RGPU=unit1 \
RGSM=$SAMPLE

SORT=$base.sort.bam
samtools sort $RG_OUTPUT > $SORT

#Remove Duplicates
METRICS=${base}_metrics.txt
Mark_OUTPUT=$base.deduped.bam
java -jar picard.jar MarkDuplicates --REMOVE_DUPLICATES true \
--INPUT $SORT --OUTPUT $Mark_OUTPUT --METRICS_FILE $METRICS

samtools index $Mark_OUTPUT

#remove intermediate files
rm $RG_OUTPUT
rm $SORT



bam=$Mark_OUTPUT
#Filter and get length per fragment
samtools view -f 3 -F 3852 -q 30 $bam | \
        awk '{if($9>0) print $9}' > ${OUTDIR}/${SAMPLE}_lengths.txt
        
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
        
    }' > ${OUTDIR}/${SAMPLE}.bed

# echo getting reads with lengths below 40bp ...
# 
# samtools view -h -f 3 -F 3852 "$bam" | \
# awk 'BEGIN {OFS="\t"} 
#      /^@/ {print; next} 
#      ($9 > 0 && $9 < 40)' >  ${OUTPUT_DIR}/${SAMPLE}_short.sam

echo "Extracting fragment lengths for $SAMPLE ..."

Rscript scripts/fragment_lengths.R "${OUTDIR}/${SAMPLE}_lengths.txt" "${OUTDIR}/${SAMPLE}_length_freq.tsv"


