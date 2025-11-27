#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor


HERE = Path(__file__).resolve().parent
print("Path is ",HERE)

def run_cmd(cmd, dry_run=False):
    """Run a shell command, print it first."""
    print("[CMD]", " ".join(str(c) for c in cmd))
    if dry_run:
        return
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(f"[ERROR] Command failed with exit code {e.returncode}", file=sys.stderr)
        sys.exit(e.returncode)


def cmd_prepare_TSS(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gtf = Path(args.gtf).resolve()
    chrom_sizes = Path(args.chrom_sizes).resolve()

    cmd = [
        str(HERE / "get_TSS.sh"),
        "--gtf", str(gtf),
        "--chrom-sizes", str(chrom_sizes),
        "--window", str(args.window),
        "--outdir", str(outdir)
    ]
    if args.genelist is not None:
        cmd += ["--genelist", str(Path(args.genelist).resolve())]

    run_cmd(cmd, dry_run=args.dry_run)

def cmd_preprocess_bam(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    bam = Path(args.bam).resolve()

    cmd = [
        str(HERE / "preprocess_bam.sh"),
        "--bam", str(bam),
        "--outdir", str(outdir),
        "--sample", args.sample
    ]
    if args.min_len is not None:
        cmd += ["--min-len", str(args.min_len)]
    if args.max_len is not None:
        cmd += ["--max-len", str(args.max_len)]

    run_cmd(cmd, dry_run=args.dry_run)
    
def cmd_coverage(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fragments = Path(args.fragments).resolve()
    regions = Path(args.regions).resolve()

    cmd = [
        str(HERE / "compute_coverage.sh"),
        "--fragments", str(fragments),
        "--regions", str(regions),
        "--outdir", str(outdir),
        "--sample", args.sample
    ]
    if args.upstream_bp is not None:
        cmd += ["--upstream-bp", str(args.upstream_bp)]
    if args.downstream_bp is not None:
        cmd += ["--downstream-bp", str(args.downstream_bp)]
    if args.len_freq is not None:
        cmd += ["--len-freq", str(args.len_freq)]
    if args.chrom_sizes is not None:
        cmd += ["--chrom-sizes", str(args.chrom_sizes)]   

    run_cmd(cmd, dry_run=args.dry_run)

def cmd_merge_plots(args):
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    cmd = [
        str(HERE / "merge_plots.sh"),
        "--sample-file", args.sample_file,
        "--outdir", str(outdir),
    ]

    run_cmd(cmd, dry_run=args.dry_run)

def cmd_batch(args):
    cmd = [
        str(HERE / "batch.sh"),
        "--sample-file", str(args.sample_file),
        "--regions", str(args.regions),
        "--upstream-bp", str(args.upstream_bp),
        "--downstream-bp", str(args.downstream_bp),
        "--chrom-sizes", str(args.chrom_sizes),

    ]

    run_cmd(cmd, dry_run=args.dry_run)



# ------------ Argument parser ------------ #

def build_parser():
    parser = argparse.ArgumentParser(
        prog="cf_frag",
        description="Command line wrapper for cfDNA fragmentomics pipeline"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="Print commands without executing them"
    )

    subparsers = parser.add_subparsers(dest="command", required=True)

    # prepare_TSS
    p_ref = subparsers.add_parser(
        "prepare_TSS",
        help="Generate TSS/region BED files from annotation and chrom sizes"
    )
    p_ref.add_argument("--gtf", required=True, help="Annotation GTF (e.g. gencode.v49.annotation.gtf)")
    p_ref.add_argument("--chrom-sizes", required=True, help="Chromosome sizes file (hg38.chrom.sizes)")
    p_ref.add_argument("--window", type=int, default=1000,
                       help="Half-window around TSS in bp (default: 2500)")
    p_ref.add_argument("--genelist", help="Optional file with gene IDs to restrict TSS set")
    p_ref.add_argument("--outdir", required=True, help="Output directory")
    p_ref.set_defaults(func=cmd_prepare_TSS)
    
    # preprocess-bam
    p_pre = subparsers.add_parser(
        "preprocess-bam",
        help="Deduplicate BAM, filter, and generate fragments BED"
    )
    p_pre.add_argument("--bam", required=True, help="Input BAM")
    p_pre.add_argument("--sample", required=True, help="Sample name (used in output filenames)")
    p_pre.add_argument("--min-len", type=int, help="Minimum fragment length")
    p_pre.add_argument("--max-len", type=int, help="Maximum fragment length")
    p_pre.add_argument("--outdir", required=True, help="Output directory")
    p_pre.set_defaults(func=cmd_preprocess_bam)
    
# coverage
    p_cov = subparsers.add_parser(
        "coverage",
        help="Compute coverage/median/normalized  over regions"
    )
    p_cov.add_argument("--fragments", required=True, help="Fragments BED (e.g. fragments.sorted.bed)")
    p_cov.add_argument("--regions", required=True, help="Regions BED (centered or TSS windows)")
    p_cov.add_argument("--sample", required=True, help="Sample name")
    p_cov.add_argument("--upstream-bp", type=int, default=2500,
                       help="Upstream flank size (bp) for normalization")
    p_cov.add_argument("--downstream-bp", type=int, default=2500,
                       help="Downstream flank size (bp) for normalization")
    p_cov.add_argument("--outdir", required=True, help="Output directory")
    p_cov.add_argument("--chrom-sizes", required=True, help="chromosome sizes")
    p_cov.add_argument("--len-freq", required=False, help="length frequncy; tsv file")
    p_cov.set_defaults(func=cmd_coverage)

# merge_plots
    p_cov = subparsers.add_parser(
        "merge-plots",
        help="Collapes frequncy histograms"
    )
    p_cov.add_argument("--sample-file", required=True, help="a text file; first column sample name second the group(e.g case,control)")
    p_cov.add_argument("--outdir", required=True, help="Output directory")
    p_cov.set_defaults(func=cmd_merge_plots)    

# batch
    # batch
    p_batch = subparsers.add_parser("batch",help="run in batch mode via batch.sh")
    p_batch.add_argument("--sample-file",required=True,help="text file; first column sample name, second the group (e.g. case, control)")
    p_batch.add_argument("--regions",required=True,help="BED file of selected regions")
    p_batch.add_argument("--upstream-bp",type=int,default=2500,help="bp upstream of regions for normalization (default 2500)")
    p_batch.add_argument("--downstream-bp",type=int,default=2500,help="bp downstream of regions for normalization (default 2500)")
    p_batch.add_argument("--chrom-sizes",required=True,help="reference genome chrom sizes")

    p_batch.set_defaults(func=cmd_batch)

    
    
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    if not hasattr(args, "func"):
        parser.print_help()
        sys.exit(1)
    args.func(args)


if __name__ == "__main__":
    main()
