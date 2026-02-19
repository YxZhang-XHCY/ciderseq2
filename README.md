# CIDER-Seq2 v3.0.0

CIDER-Seq2 is a computational pipeline for detecting **extrachromosomal circular DNA (eccDNA)** from PacBio HiFi long-read sequencing data. It takes advantage of the **Rolling Circle Amplification (RCA)** signature in HiFi reads — circular DNA molecules are amplified into tandem repeat concatemers during SMRT library preparation, and CIDER-Seq2 identifies these tandem repeats to recover the original circular sequences and map them to a host reference genome.

**v3.0.0** replaces the original MUSCLE-based DeConcat algorithm with [TideHunter](https://github.com/yangao07/TideHunter) for faster and more accurate RCA tandem repeat detection, and packages the entire pipeline as a pip-installable Python package with a single `ciderseq2` command.

## Table of Contents

- [How It Works](#how-it-works)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Usage](#usage)
  - [Full Pipeline (Single Sample)](#full-pipeline-single-sample)
  - [Multiple Biological Replicates](#multiple-biological-replicates)
  - [Skip DeConcat Mode](#skip-deconcat-mode)
  - [Standalone Sub-commands](#standalone-sub-commands)
- [Parameters](#parameters)
- [Output Files](#output-files)
- [SLURM Batch Job Example](#slurm-batch-job-example)
- [Pipeline Details](#pipeline-details)
- [Dependencies](#dependencies)
- [FAQ](#faq)
- [Citation](#citation)
- [License](#license)

## How It Works

```
HiFi reads (FASTA)          Reference genome (FASTA)
        │                            │
        ▼                            │
┌──────────────────┐                 │
│ Step 1: DeConcat │                 │
│   (TideHunter)   │                 │
└────────┬─────────┘                 │
         │ consensus sequences       │
         │ + copy number stats       │
         ▼                           ▼
┌──────────────────────────────────────┐
│ Step 2: eccDNA Detection             │
│   - Filter reads with copy > 0      │
│   - BLAST against host genome       │
│   - Reorder partial BLAST hits      │
│   - Assess single-hit coverage      │
│   - Cross-replicate filtering (opt) │
└────────┬─────────────────────────────┘
         │
         ▼
┌──────────────────┐
│ Step 3: BED      │
│   Generation     │
└────────┬─────────┘
         │
         ▼
  eccDNA candidates (FASTA)
  Genomic coordinates (BED)
```

**Step 1 — DeConcat (TideHunter):** Each HiFi read is scanned for tandem repeats. If a read contains >= 2 copies of a repeated unit (indicating RCA origin), TideHunter extracts the consensus sequence. Non-circular reads pass through unchanged. The output includes a consensus FASTA and a statistics file recording the copy number per read.

**Step 2 — eccDNA Detection:** Reads with copy number > 0 are BLASTed against the host reference genome. Partial BLAST hits are merged (reordered) to handle cases where a circular sequence spans two adjacent genomic hits. Reads whose reordered sequence maps to a single genomic locus with >= 95% query coverage are classified as eccDNA candidates. For multi-replicate experiments, candidates are further filtered to retain only those found in multiple samples.

**Step 3 — BED Generation:** The genomic coordinates of each eccDNA candidate are extracted from the BLAST results and written as a BED file with length and copy number annotations.

## Installation

### 1. Create conda environment

```bash
conda create -n ciderseq2 python=3.7
conda activate ciderseq2
```

### 2. Install external dependencies

```bash
conda install -c bioconda tidehunter blast
pip install biopython==1.77 click
```

> **Note:** Biopython must be version 1.77 or earlier. The `Bio.Blast.Applications` module used by the eccDNA detection step was removed in Biopython 1.78+.

### 3. Install ciderseq2

**Option A — Local install (recommended for servers without GitHub access)**

```bash
# On a machine with GitHub access:
git clone https://github.com/YxZhang-XHCY/ciderseq2.git
cd ciderseq2
git checkout tidehunter-deconcat

# Transfer to server:
scp -r ciderseq2/ user@server:/path/to/tools/

# On the server:
conda activate ciderseq2
cd /path/to/tools/ciderseq2
pip install .
```

**Option B — Direct install from GitHub**

```bash
pip install git+https://github.com/YxZhang-XHCY/ciderseq2.git@tidehunter-deconcat
```

### 4. Verify installation

```bash
ciderseq2 --version    # 3.0.0
ciderseq2 --help
ciderseq2-deconcat --help
ciderseq2-to-bed --help
```

## Quick Start

```bash
conda activate ciderseq2

ciderseq2 \
    --input sample.hifi.fasta \
    --genome reference_genome.fasta \
    --output results/ \
    --prefix my_run \
    --threads 8
```

This will:
1. Run TideHunter to identify circular reads and extract consensus sequences
2. BLAST consensus sequences against the reference genome to detect eccDNA
3. Generate a BED file with genomic coordinates of each eccDNA

## Usage

### Full Pipeline (Single Sample)

```bash
ciderseq2 \
    --input sample.hifi.fasta \
    --genome reference_genome.fasta \
    --output results/ \
    --prefix my_run \
    --threads 8
```

### Multiple Biological Replicates

Provide comma-separated input files. The pipeline will process each sample independently, then perform cross-replicate filtering to identify shared eccDNA:

```bash
ciderseq2 \
    --input rep1.fasta,rep2.fasta,rep3.fasta \
    --genome genome.fasta \
    --output results/ \
    --prefix multi_rep \
    --threads 8
```

### Skip DeConcat Mode

If you already have deconcat output from a previous run or from another tool, skip Step 1 with `--skip-deconcat`. The input format changes to `fasta:stat` pairs:

```bash
# Single sample
ciderseq2 \
    --input deconcat.fa:deconcat.stat \
    --genome genome.fasta \
    --output results/ \
    --prefix rerun \
    --skip-deconcat

# Multiple replicates
ciderseq2 \
    --input rep1.fa:rep1.stat,rep2.fa:rep2.stat \
    --genome genome.fasta \
    --output results/ \
    --prefix rerun_multi \
    --skip-deconcat
```

### Standalone Sub-commands

Run individual steps independently:

```bash
# Step 1 only: TideHunter deconcat
ciderseq2-deconcat \
    --input raw.hifi.fasta \
    --output-fa out.deconcat.fa \
    --output-stat out.deconcat.stat \
    --threads 8

# Step 3 only: Generate BED from existing eccDNA results
ciderseq2-to-bed \
    --assess output/assess1.fasta \
    --blast output/reorder1.fasta_on_genomedb.blastn6 \
    --stat deconcat.stat \
    --output eccdna.bed
```

## Parameters

### `ciderseq2` (full pipeline)

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | HiFi FASTA file(s). Comma-separated for multiple replicates. With `--skip-deconcat`, use `fasta:stat` format |
| `--genome` | *required* | Host reference genome FASTA |
| `--output` | *required* | Output root directory |
| `--prefix` | *required* | Output prefix (creates `<output>/<prefix>/` subdirectory) |
| `--threads` | 8 | Number of threads for TideHunter and BLAST |
| `--gap-window` | 150 | Gap window (bp) for merging adjacent partial BLAST hits |
| `--tidehunter` | auto-detect | Path to TideHunter executable |
| `--skip-deconcat` | off | Skip Step 1; input must be pre-computed `fasta:stat` pairs |

### `ciderseq2-deconcat`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | *required* | Input HiFi FASTA |
| `--output-fa` | *required* | Output deconcat FASTA (consensus for circular reads, original for linear) |
| `--output-stat` | *required* | Output statistics file (11-column TSV) |
| `--tidehunter` | auto-detect | Path to TideHunter executable |
| `--threads` | 8 | Number of threads |
| `--min-copy` | 2 | Minimum tandem repeat copy number |
| `--min-period` | 30 | Minimum consensus length (bp) |
| `--max-period` | 10000 | Maximum consensus length (bp) |

### `ciderseq2-to-bed`

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--assess` | *required* | Assessed eccDNA FASTA file(s) (repeatable: `--assess f1 --assess f2`) |
| `--blast` | *required* | Reorder BLAST6 result(s) (matching `--assess` order) |
| `--output` | *required* | Output BED file |
| `--stat` | optional | Deconcat stat file(s) for copy number annotation |

## Output Files

### Directory structure

```
<output>/<prefix>/
├── 01_deconcat/
│   ├── <sample>.deconcat.fa            # Consensus sequences (circular reads replaced)
│   └── <sample>.deconcat.stat          # Per-read statistics (11-column TSV)
│
├── 02_eccDNA/
│   ├── samples.list                    # Input list for eccDNA detection
│   └── output/
│       ├── <prefix>_eccDNA_candidates.fasta   ★ Final eccDNA sequences
│       ├── <prefix>_eccDNA.bed                ★ Genomic coordinates
│       ├── assess*.fasta                  # Per-sample eccDNA candidates
│       ├── reorder*.fasta                 # Reordered reads (partial hits merged)
│       ├── reorder*.fasta_on_genomedb.blastn6  # BLAST results
│       ├── rfile*.fasta                   # Filtered reads (copy number > 0)
│       ├── replicatefile*.fasta           # Cross-replicate shared eccDNA (multi-sample only)
│       └── genomedb.*                     # BLAST database files
│
└── logs/                                  # Pipeline logs
```

### Deconcat stat file format (11-column TSV)

| Column | Description |
|--------|-------------|
| 1 | Read ID |
| 2 | Copy number (0 = linear, >= 2 = circular) |
| 3 | Average match score |
| 4-11 | Reserved (set to 0) |

### BED file format (8 columns)

| Column | Description |
|--------|-------------|
| chrom | Chromosome name |
| start | Start position (0-based) |
| end | End position |
| name | Read ID |
| score | Alignment identity x 10 (0-1000) |
| strand | + or - |
| length | eccDNA length in bp |
| copies | RCA tandem repeat copy number |

## SLURM Batch Job Example

For running multiple samples on a computing cluster:

```bash
#!/bin/bash
#SBATCH --job-name=ciderseq2
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err
#SBATCH --cpus-per-task=16
#SBATCH --mem=100G
#SBATCH --time=96:00:00
#SBATCH --array=0-2

eval "$(conda shell.bash hook)"
conda activate ciderseq2

SAMPLES=("sample1.hifi.fasta" "sample2.hifi.fasta" "sample3.hifi.fasta")
GENOME="/path/to/reference_genome.fasta"
OUTPUT="/path/to/results"

SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}
PREFIX=$(basename ${SAMPLE%.hifi.fasta})

ciderseq2 \
    --input ${SAMPLE} \
    --genome ${GENOME} \
    --output ${OUTPUT} \
    --prefix ${PREFIX} \
    --threads ${SLURM_CPUS_PER_TASK}
```

Submit with:
```bash
mkdir -p logs
sbatch run_ciderseq2.sh
```

## Pipeline Details

### Step 1: DeConcat (TideHunter)

[TideHunter](https://github.com/yangao07/TideHunter) scans each HiFi read for tandem repeats using a partial-order alignment approach. For each read:

- If **>= 2 tandem repeat copies** are found: the read is classified as circular. TideHunter extracts the consensus sequence of the repeated unit, which represents the original eccDNA sequence. The copy number is recorded.
- If **no tandem repeat** is found: the read passes through unchanged with copy number = 0.

This replaces the original MUSCLE-based DeConcat algorithm from CIDER-Seq v1/v2, providing significantly faster processing for long HiFi reads (typically 10-20 kb).

### Step 2: eccDNA Detection

The eccDNA detection module (derived from `cs-eccDNA.py` by Cornet & Zaidi) performs:

1. **Filtering:** Only reads with copy number > 0 (circular) are retained.
2. **BLAST:** Filtered reads are BLASTed against the host reference genome.
3. **Reordering:** Partial BLAST hits on the same chromosome within the gap window are merged. This handles cases where a circular sequence is linearized at a position that splits a genomic alignment into two adjacent hits.
4. **Assessment:** Reordered reads are re-BLASTed. If the best hit covers >= 95% of the query length, the read is classified as an eccDNA candidate.
5. **Cross-replicate filtering** (multi-sample only): eccDNA candidates found in multiple biological replicates are identified by all-vs-all BLAST at >= 99% identity.

### Step 3: BED Generation

For each eccDNA candidate, the primary BLAST hit (best bitscore) determines the genomic coordinates. The output BED file includes length and copy number annotations for downstream analysis.

## Dependencies

| Dependency | Version | Purpose | Install |
|------------|---------|---------|---------|
| Python | >= 3.7 | Runtime | conda |
| [TideHunter](https://github.com/yangao07/TideHunter) | >= 1.5 | Tandem repeat detection (Step 1) | `conda install -c bioconda tidehunter` |
| [BLAST+](https://blast.ncbi.nlm.nih.gov/) | >= 2.10 | Genome alignment (Step 2) | `conda install -c bioconda blast` |
| [Biopython](https://biopython.org/) | <= 1.77 | BLAST command wrapper | `pip install biopython==1.77` |
| [Click](https://click.palletsprojects.com/) | any | CLI framework | `pip install click` |

All Python dependencies (`biopython`, `click`) are automatically installed by `pip install .`

External tools (`TideHunter`, `blastn`, `makeblastdb`) must be installed separately via conda/bioconda and available in `$PATH`.

## FAQ

**Q: Why does Biopython need to be <= 1.77?**

The eccDNA detection module uses `Bio.Blast.Applications.NcbiblastnCommandline`, which was removed in Biopython 1.78. Use `pip install biopython==1.77` to ensure compatibility.

**Q: Can I use nanopore reads instead of HiFi?**

The pipeline assumes PacBio HiFi reads with RCA tandem repeat structure. Nanopore reads from RCA-based protocols may work, but this has not been tested. The key requirement is that circular DNA reads contain >= 2 tandem copies of the original sequence.

**Q: What is the `--gap-window` parameter?**

When a circular sequence is linearized at a breakpoint that falls within a genomic region, BLAST may report two adjacent hits instead of one. The gap window (default: 150 bp) controls the maximum allowed gap between two hits on the same chromosome to merge them into a single eccDNA candidate.

**Q: How does multi-replicate filtering work?**

With >= 2 input samples, after per-sample eccDNA identification, candidates are compared across samples using all-vs-all BLAST. Reads matching at >= 99% identity between replicates are retained. This reduces false positives from random linear DNA fragments.

**Q: My server can't access GitHub. How do I install?**

Clone the repository on a machine with internet access, then transfer:
```bash
# Local machine
git clone https://github.com/YxZhang-XHCY/ciderseq2.git
cd ciderseq2 && git checkout tidehunter-deconcat

# Transfer to server
scp -r ciderseq2/ user@server:/path/to/tools/

# On server
cd /path/to/tools/ciderseq2 && pip install .
```

## Citation

Mehta D, Cornet L, Hirsch-Hoffmann M, Zaidi SSA, Vanderschuren H (2020) _Full-length sequencing of circular DNA viruses and extra-chromosomal circular DNA using CIDER-Seq._ **Nature Protocols**, 15, 1673-1689. [doi:10.1038/s41596-020-0301-0](https://doi.org/10.1038/s41596-020-0301-0)

Mehta D, Hirsch-Hoffmann M, Were M, Patrignani A, Were H, Gruissem W, Vanderschuren H (2019) _A new full-length circular DNA sequencing method for viral-sized genomes reveals that RNAi transgenic plants provoke a shift in geminivirus populations in the field._ **Nucleic Acids Research**, 47(2). [doi:10.1093/nar/gky914](https://doi.org/10.1093/nar/gky914)

## License

Copyright &copy; 2019, Matthias Hirsch-Hoffmann and Devang Mehta. Licensed under the [GNU General Public License v3.0](LICENSE.txt).
