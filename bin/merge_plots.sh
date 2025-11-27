#!/bin/bash

# Parse command-line arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --sample-file) sample_file="$2"; shift ;;
        --outdir) OUTDIR="$2"; shift ;;
        *) 
            echo "Unknown parameter passed: $1"; 
            exit 1 ;;
    esac
    shift
done

Rscript scripts/merge.R $sample_file $OUTDIR