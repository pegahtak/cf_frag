# cfDNA Fragmentomics Pipeline

Command-line pipeline for basic cfDNA fragmentomics analysis, including:

- extracting TSS-centered regions
- preprocessing BAM files into fragment BEDs
- computing coverage and length-frequency profiles
- generating per-sample plots and merged case/control histograms
- running all steps in batch for multiple samples

All steps are wrapped by a single CLI: `bin/cf_frags.py`.











This pipeline inputs bam files, sorts them and produces fragmentomics 
features and plot them.


Test data obtained from cfDNA pipe paper. both test and control

Optional
#prepare Transcription Start Sites.
gene annotations of GRCh38 genome has been obtained from GENCODE and placed in the data folder.(https://www.gencodegenes.org/human/)
You can replace other versions with --gtf option.
To obtain the TSS regions run 
cd bin
cf_frags.py prepare_TSS --gtf data/gencode.v49.annotation.gtf \
--chrom-sizes data/hg38.chrom.sizes --outdir data
This code obtains start point of every transcript and then collapeses TSS for transcripts with the same gene id
then  a window of +/- 1kb around the start point is considered as the TSS.
You can change the size of window by the option --window (e.g. --window 2500)

#prepare BAM files


