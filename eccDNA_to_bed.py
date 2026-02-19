#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Generate BED file from CIDER-Seq2 eccDNA detection results.

Parses the reorder BLAST results and filters by assessed eccDNA read IDs
to produce a BED6 file with genomic coordinates of each eccDNA.

For each assessed read, takes the primary BLAST hit (first hit per read,
which is the best by bitscore) to determine the genomic origin.
"""

import sys
import click


def parse_assess_ids(assess_fasta):
    """Extract read IDs from assess FASTA."""
    ids = set()
    with open(assess_fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.add(line.strip()[1:].split()[0])
    return ids


def parse_blast6_primary(blast6_path, valid_ids):
    """Parse blastn6, keep only primary hit (first per query) for valid IDs.

    Returns list of (chrom, start, end, name, score, strand).
    """
    records = []
    seen = set()

    with open(blast6_path) as f:
        for line in f:
            cols = line.strip().split("\t")
            if len(cols) < 12:
                continue
            query = cols[0]
            if query in seen:
                continue
            if query not in valid_ids:
                continue
            seen.add(query)

            chrom = cols[1]
            identity = float(cols[2])
            aln_len = int(cols[3])
            sstart = int(cols[8])
            send = int(cols[9])

            # BED is 0-based half-open
            if sstart <= send:
                bed_start = sstart - 1
                bed_end = send
                strand = "+"
            else:
                bed_start = send - 1
                bed_end = sstart
                strand = "-"

            score = int(identity * 10)  # scale identity to 0-1000
            records.append((chrom, bed_start, bed_end, query, score, strand))

    return records


@click.command()
@click.option("--assess", required=True, multiple=True,
              help="Assessed eccDNA FASTA (can specify multiple)")
@click.option("--blast", required=True, multiple=True,
              help="Reorder BLAST6 result (matching --assess order)")
@click.option("--output", required=True, help="Output BED file")
@click.option("--stat", multiple=True, default=None,
              help="Deconcat stat file(s) for copy number annotation")
def main(assess, blast, output, stat):
    """Generate BED file from CIDER-Seq2 eccDNA results."""

    if len(assess) != len(blast):
        sys.exit("Number of --assess and --blast files must match")

    # Load copy number from stat files if provided
    copy_of = {}
    for stat_file in stat:
        with open(stat_file) as f:
            for line in f:
                cols = line.strip().split("\t")
                if len(cols) >= 2 and int(cols[1]) > 0:
                    copy_of[cols[0]] = cols[1]

    all_records = []
    for assess_fa, blast6 in zip(assess, blast):
        ids = parse_assess_ids(assess_fa)
        records = parse_blast6_primary(blast6, ids)
        all_records.extend(records)

    # Sort by chrom, start
    all_records.sort(key=lambda r: (r[0], r[1]))

    with open(output, "w") as out:
        for chrom, start, end, name, score, strand in all_records:
            length = end - start
            copies = copy_of.get(name, ".")
            # BED6 + extra columns: length, copy_number
            out.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}"
                      f"\t{length}\t{copies}\n")

    click.echo(f"Wrote {len(all_records)} eccDNA records to {output}", err=True)


if __name__ == "__main__":
    main()
