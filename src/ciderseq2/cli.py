#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
CIDER-Seq2 eccDNA Detection Pipeline — main CLI entry point.

Replaces ciderseq2_pipeline.sh with a pure-Python Click CLI.
"""

import glob
import os
import shutil
import subprocess
import sys

import click

from ciderseq2 import __version__
from ciderseq2.deconcat import run_deconcat
from ciderseq2.to_bed import generate_bed
from ciderseq2.eccDNA.pipeline import run_eccdna


def check_executable(name):
    """Check if an external command is available in PATH."""
    path = shutil.which(name)
    if path is None:
        click.echo(f"ERROR: Required command not found: {name}", err=True)
        click.echo("Please install via: conda install -c bioconda tidehunter blast", err=True)
        sys.exit(1)
    return path


@click.command()
@click.option("--input", "input_files", required=True,
              help="Input HiFi FASTA file(s), comma-separated for replicates")
@click.option("--genome", required=True, type=click.Path(exists=True),
              help="Reference genome FASTA")
@click.option("--output", "outdir", required=True,
              help="Output directory")
@click.option("--prefix", required=True,
              help="Output prefix")
@click.option("--threads", default=8, show_default=True,
              help="Number of BLAST threads")
@click.option("--gap-window", default=150, show_default=True,
              help="Gap window for merging partial BLAST hits")
@click.option("--tidehunter", "tidehunter_exe", default=None,
              help="Path to TideHunter executable (default: auto-detect)")
@click.option("--skip-deconcat", is_flag=True, default=False,
              help="Skip DeConcat; input format is fasta:stat[,fasta:stat,...]")
@click.version_option(__version__)
def main(input_files, genome, outdir, prefix, threads, gap_window,
         tidehunter_exe, skip_deconcat):
    """CIDER-Seq2: eccDNA detection from HiFi long reads."""

    genome = os.path.abspath(genome)

    # ── Check dependencies ────────────────────────────────────────────
    click.echo(f"CIDER-Seq2 v{__version__}", err=True)
    click.echo("Checking dependencies...", err=True)

    check_executable("blastn")
    check_executable("makeblastdb")

    if tidehunter_exe is None and not skip_deconcat:
        tidehunter_exe = shutil.which("TideHunter")
        if tidehunter_exe is None:
            click.echo("ERROR: TideHunter not found in PATH. "
                        "Use --tidehunter to specify path.", err=True)
            sys.exit(1)

    # Check Biopython
    try:
        from Bio.Blast.Applications import NcbiblastnCommandline  # noqa: F401
    except ImportError:
        click.echo("ERROR: Biopython (<=1.77) required. "
                    "Install: pip install 'biopython<=1.77'", err=True)
        sys.exit(1)

    click.echo("Dependencies OK", err=True)

    # ── Setup directory structure ─────────────────────────────────────
    workdir = os.path.join(outdir, prefix)
    dir_deconcat = os.path.join(workdir, "01_deconcat")
    dir_eccdna = os.path.join(workdir, "02_eccDNA")
    dir_log = os.path.join(workdir, "logs")

    os.makedirs(dir_deconcat, exist_ok=True)
    os.makedirs(os.path.join(dir_eccdna, "output"), exist_ok=True)
    os.makedirs(dir_log, exist_ok=True)

    # ── Parse input samples ───────────────────────────────────────────
    samples = input_files.split(",")
    num_samples = len(samples)
    click.echo(f"Samples: {num_samples}", err=True)

    # ── Step 1: DeConcat ──────────────────────────────────────────────
    deconcat_fastas = []
    deconcat_stats = []

    if skip_deconcat:
        click.echo("[Step 1] Skipping DeConcat (using pre-computed outputs)", err=True)
        for sample in samples:
            if ":" not in sample:
                click.echo(f"ERROR: --skip-deconcat requires fasta:stat format, got: {sample}",
                            err=True)
                sys.exit(1)
            fasta, stat = sample.split(":", 1)
            fasta = os.path.abspath(fasta)
            stat = os.path.abspath(stat)
            if not os.path.isfile(fasta):
                sys.exit(f"Fasta file not found: {fasta}")
            if not os.path.isfile(stat):
                sys.exit(f"Stat file not found: {stat}")
            deconcat_fastas.append(fasta)
            deconcat_stats.append(stat)
            click.echo(f"  {fasta} | {stat}", err=True)
    else:
        click.echo("[Step 1] DeConcat (TideHunter)", err=True)
        for idx, sample in enumerate(samples, 1):
            sample_file = os.path.abspath(sample)
            if not os.path.isfile(sample_file):
                sys.exit(f"Input file not found: {sample_file}")

            sample_name = os.path.splitext(os.path.basename(sample_file))[0]
            click.echo(f"  Processing sample {idx}/{num_samples}: {sample_name}", err=True)

            clean_fa = os.path.join(dir_deconcat, f"{sample_name}.deconcat.fa")
            clean_stat = os.path.join(dir_deconcat, f"{sample_name}.deconcat.stat")

            run_deconcat(
                input_fasta=sample_file,
                output_fa=clean_fa,
                output_stat=clean_stat,
                tidehunter_exe=tidehunter_exe,
                threads=threads,
            )

            deconcat_fastas.append(clean_fa)
            deconcat_stats.append(clean_stat)

    click.echo(f"DeConcat outputs: {len(deconcat_fastas)} samples", err=True)

    # ── Step 2: eccDNA Detection ──────────────────────────────────────
    click.echo("[Step 2] eccDNA Detection", err=True)

    # Create list file (tab-separated: fasta<TAB>stat)
    list_file = os.path.join(dir_eccdna, "samples.list")
    with open(list_file, "w") as f:
        for fa, st in zip(deconcat_fastas, deconcat_stats):
            f.write(f"{fa}\t{st}\n")

    click.echo(f"  List file: {list_file}", err=True)

    # os.chdir() is required: eccDNA modules use hardcoded 'output/' paths
    saved_cwd = os.getcwd()
    os.chdir(dir_eccdna)

    try:
        run_eccdna(
            list_file=list_file,
            genome_file=genome,
            blast_threads=threads,
            gap_window=gap_window,
        )
    finally:
        os.chdir(saved_cwd)

    # ── Step 3: Summary & BED ─────────────────────────────────────────
    click.echo("[Step 3] Summarizing results", err=True)

    eccdna_output = os.path.join(dir_eccdna, "output")

    # Per-sample summary
    for i in range(1, num_samples + 1):
        assess = os.path.join(eccdna_output, f"assess{i}.fasta")
        if os.path.isfile(assess):
            n = sum(1 for line in open(assess) if line.startswith(">"))
            click.echo(f"  Sample {i} assessed eccDNA: {n}", err=True)

    # Single sample: copy assess1 as final output
    if num_samples == 1:
        assess1 = os.path.join(eccdna_output, "assess1.fasta")
        if os.path.isfile(assess1):
            final_out = os.path.join(eccdna_output, f"{prefix}_eccDNA_candidates.fasta")
            shutil.copy2(assess1, final_out)
            click.echo(f"  Final candidates: {final_out}", err=True)

    # Generate BED file
    click.echo("Generating BED file...", err=True)

    bed_assess = []
    bed_blast = []
    for i in range(1, num_samples + 1):
        assess = os.path.join(eccdna_output, f"assess{i}.fasta")
        blast = os.path.join(eccdna_output, f"reorder{i}.fasta_on_genomedb.blastn6")
        if os.path.isfile(assess) and os.path.isfile(blast):
            bed_assess.append(assess)
            bed_blast.append(blast)

    if bed_assess:
        bed_out = os.path.join(eccdna_output, f"{prefix}_eccDNA.bed")
        generate_bed(
            assess_files=bed_assess,
            blast_files=bed_blast,
            output=bed_out,
            stat_files=deconcat_stats,
        )

    click.echo("=" * 60, err=True)
    click.echo(f"Results: {eccdna_output}/", err=True)
    click.echo(f"Logs:    {dir_log}/", err=True)
    click.echo("Pipeline complete!", err=True)


if __name__ == "__main__":
    main()
