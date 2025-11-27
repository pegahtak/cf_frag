# cfDNA Fragmentomics Pipeline

Command-line pipeline for basic cfDNA fragmentomics analysis, including:

- extracting TSS-centered regions
- preprocessing BAM files into fragment BEDs
- computing coverage and length-frequency profiles
- generating per-sample plots and merged case/control histograms
- running all steps in batch for multiple samples

All steps are wrapped by a single CLI: `bin/cf_frags.py`.

## 1. Installation and requirements

### Dependencies

- Unix-like system (bash)
- Python ≥ 3.6
- R ≥ 4.0
- Common bioinformatics tools (used inside the bash scripts):
  - `samtools`
  - `bedtools`
  - `awk`, `gawk`,`sed`, `grep`, `sort`, `gzip`, etc.

R packages (used in the plotting scripts):

- `ggplot2`
-  `data.table`

### Setup

Clone this repository and make the scripts executable:

```bash
https://github.com/pegahtak/cf_frag.git
cd fragmentomics
```
make the scripts executable:
```bash
chmod +x bin/cf_frags.py
chmod +x bin/get_TSS.sh
chmod +x bin/preprocess_bam.sh
chmod +x bin/compute_coverage.sh
chmod +x bin/merge_plots.sh
chmod +x bin/batch.sh
```

## 2. Directory layout (expected)

A typical layout:

```text
fragmentomics/
├── bin/
│   ├── cf_frags.py          # main CLI wrapper
│   ├── get_TSS.sh
│   ├── preprocess_bam.sh
│   ├── compute_coverage.sh
│   └── merge_plots.sh
├── scripts/
│   └── merge.R              # R script that draws merged histograms
├── data/
│   ├── gencode.vXX.annotation.gtf
│   ├── hg38.chrom.sizes
│   ├── manifest.tsv         # sample → BAM mapping
│   └── sample_files.txt     # sample → group mapping
└── results/
    └── ...                  # created by the pipeline
```







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


