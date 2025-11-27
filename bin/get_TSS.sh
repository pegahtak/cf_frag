#!/bin/bash

# Parse command-line arguments
GTF=""
CHRM_SIZES=""
WINDOW=1000
OUTDIR=""
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --gtf) GTF="$2"; shift ;;
        --chrom-sizes) CHRM_SIZES="$2"; shift ;;
        --window) WINDOW="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

# Check required arguments
if [[ -z "$GTF" || -z "$CHRM_SIZES" || -z "$OUTDIR" ]]; then
    echo "--gtf  --chrom-sizes  --outdir are required"
    exit 1
fi

awk 'BEGIN{FS=OFS="\t"}
$3=="transcript"{
  gene_id="";
  transcript_id="";

  # parse field 9: key "value"; key "value";
  n = split($9, a, ";");
  for (i=1; i<=n; i++){
    gsub(/^ +/,"",a[i]);      # trim leading spaces
    if (a[i] ~ /^gene_id /){
      split(a[i], b, "\"");
      gene_id = b[2];
    } else if (a[i] ~ /^transcript_id /){
      split(a[i], b, "\"");
      transcript_id = b[2];
    }
  }

  if (gene_id=="" || transcript_id=="") next;

  if ($7=="+"){
    tss = $4;   # 1-based
  } else if ($7=="-"){
    tss = $5;   # 1-based
  } else {
    next;
  }

  # store TSS as 1-based in both start and end columns for now
  print $1, tss, tss, gene_id ";" transcript_id, 0, $7;
}' $GTF > $OUTDIR/tss_transcripts.tmp


awk 'BEGIN{FS=OFS="\t"}
{
  split($4, ids, ";");
  gene   = ids[1];
  strand = $6;
  tss    = $2;   # 1-based TSS
  chr    = $1;

  key = gene FS strand;

  if (!(key in best_tss)){
    best_tss[key] = tss;
    best_chr[key] = chr;
  } else {
    cur = best_tss[key];
    if (strand=="+" && tss < cur){
      best_tss[key] = tss;
      best_chr[key] = chr;
    } else if (strand=="-" && tss > cur){
      best_tss[key] = tss;
      best_chr[key] = chr;
    }
  }
}
END{
  for (key in best_tss){
    split(key, ks, FS);
    gene   = ks[1];
    strand = ks[2];
    tss1   = best_tss[key];   # 1-based
    chr    = best_chr[key];

    # convert to 0-based [start,end) point
    tss0_start = tss1 - 1;
    tss0_end   = tss1;

    if (tss0_start < 0) tss0_start = 0;

    print chr, tss0_start, tss0_end, gene, 0, strand;
  }
}' $OUTDIR/tss_transcripts.tmp > $OUTDIR/gene_TSS_point.bed


#extend TSS WINDOW bp 
#check if windows are in chromosome intervals
awk -v WINDOW="$WINDOW" 'BEGIN{FS=OFS="\t"}
NR==FNR {
  # first file: chrom sizes
  chrLen[$1] = $2
  next
}
{
  chr    = $1
  tss0   = $2   # 0-based
  gene   = $4
  strand = $6

  start = tss0 - WINDOW
  end   = tss0 + WINDOW

  if (start < 0) start = 0

  # clamp to chromosome length (BED end is exclusive)
  if (chr in chrLen) {
    if (end > chrLen[chr]) end = chrLen[chr]
  }

  print chr, start, end, gene, 0, strand
}' $CHRM_SIZES $OUTDIR/gene_TSS_point.bed > $OUTDIR/gene_TSS_window.bed

#remove tmep file
rm $OUTDIR/tss_transcripts.tmp