#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
TideHunter-based DeConcat replacement for CIDER-Seq2.

Replaces the MUSCLE-based DeConcat step with TideHunter for
detecting RCA tandem repeats in HiFi long reads.

Output formats are compatible with cs-eccDNA.py (Step 2).
"""

import subprocess
import sys
import shutil
from collections import OrderedDict

import click


def parse_fasta(fasta_path):
    """Parse FASTA file into OrderedDict {header: sequence}."""
    sequences = OrderedDict()
    current_header = None
    current_seq = []

    with open(fasta_path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header is not None:
                    sequences[current_header] = "".join(current_seq)
                current_header = line[1:].split()[0]  # first word only, no spaces
                current_seq = []
            else:
                current_seq.append(line)
        if current_header is not None:
            sequences[current_header] = "".join(current_seq)

    return sequences


def run_tidehunter(input_fasta, tidehunter_exe, threads, min_copy, min_period, max_period):
    """Run TideHunter and return parsed results.

    Returns dict: {readName: (copyNum, aveMatch, consSeq)}
    """
    cmd = [
        tidehunter_exe,
        "-f", "2",       # tabular output
        "-l",            # longest tandem repeat only
        "-t", str(threads),
        "-c", str(min_copy),
        "-p", str(min_period),
        "-P", str(max_period),
        input_fasta,
    ]

    click.echo(f"Running: {' '.join(cmd)}", err=True)

    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        click.echo(f"TideHunter stderr:\n{result.stderr}", err=True)
        sys.exit(f"TideHunter failed with exit code {result.returncode}")

    # Parse tabular output (11 columns):
    # readName repN copyNum readLen start end consLen aveMatch fullLen subPos consSeq
    hits = {}
    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        cols = line.split("\t")
        if len(cols) < 11:
            continue
        read_name = cols[0]
        copy_num = float(cols[2])
        ave_match = float(cols[7])
        cons_seq = cols[10]
        hits[read_name] = (copy_num, ave_match, cons_seq)

    return hits


def run_deconcat(input_fasta, output_fa, output_stat, tidehunter_exe=None,
                 threads=8, min_copy=2, min_period=30, max_period=10000):
    """Run DeConcat as a callable function (no Click context needed).

    Args:
        input_fasta: Path to input FASTA file.
        output_fa: Path to output deconcat FASTA.
        output_stat: Path to output deconcat stat TSV.
        tidehunter_exe: Path to TideHunter executable (auto-detect if None).
        threads: Number of threads.
        min_copy: Minimum copy number for tandem repeat.
        min_period: Minimum period size.
        max_period: Maximum period size.
    """
    # Auto-detect TideHunter
    if tidehunter_exe is None:
        tidehunter_exe = shutil.which("TideHunter")
        if tidehunter_exe is None:
            sys.exit("TideHunter not found in PATH. Use --tidehunter to specify path.")
    click.echo(f"TideHunter: {tidehunter_exe}", err=True)

    # Parse input FASTA
    click.echo(f"Reading input: {input_fasta}", err=True)
    sequences = parse_fasta(input_fasta)
    click.echo(f"  Total reads: {len(sequences)}", err=True)

    # Run TideHunter
    hits = run_tidehunter(input_fasta, tidehunter_exe, threads,
                          min_copy, min_period, max_period)
    click.echo(f"  Tandem repeat reads: {len(hits)}", err=True)

    # Write outputs
    n_circular = 0
    with open(output_fa, "w") as fa_out, open(output_stat, "w") as stat_out:
        for read_name, seq in sequences.items():
            if read_name in hits:
                copy_num, ave_match, cons_seq = hits[read_name]
                rounds = int(copy_num)
                score = ave_match
                # Write consensus sequence
                fa_out.write(f">{read_name}\n{cons_seq}\n")
                n_circular += 1
            else:
                rounds = 0
                score = 0
                # Write original sequence
                fa_out.write(f">{read_name}\n{seq}\n")

            # 11-column TSV stat line
            stat_out.write(f"{read_name}\t{rounds}\t{score}\t0\t0\t0\t0\t0\t0\t0\t0\n")

    click.echo(f"Output FASTA: {output_fa} ({len(sequences)} reads)", err=True)
    click.echo(f"Output stat:  {output_stat} ({n_circular} circular)", err=True)


@click.command()
@click.option("--input", "input_fasta", required=True, help="Input FASTA file")
@click.option("--output-fa", required=True, help="Output deconcat FASTA file")
@click.option("--output-stat", required=True, help="Output deconcat stat file (11-col TSV)")
@click.option("--tidehunter", "tidehunter_exe", default=None,
              help="Path to TideHunter executable (default: auto-detect)")
@click.option("--threads", default=8, help="Number of threads")
@click.option("--min-copy", default=2, help="Minimum copy number for tandem repeat")
@click.option("--min-period", default=30, help="Minimum period size (consensus length)")
@click.option("--max-period", default=10000, help="Maximum period size (consensus length)")
def main(input_fasta, output_fa, output_stat, tidehunter_exe, threads,
         min_copy, min_period, max_period):
    """TideHunter-based DeConcat for CIDER-Seq2."""
    run_deconcat(input_fasta, output_fa, output_stat, tidehunter_exe,
                 threads, min_copy, min_period, max_period)


if __name__ == "__main__":
    main()
