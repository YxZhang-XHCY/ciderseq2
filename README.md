# CIDER-Seq2 v3.0.0

eccDNA detection from PacBio HiFi long reads.

**v3.0.0 changes:** Replaced MUSCLE-based DeConcat with [TideHunter](https://github.com/yangao07/TideHunter) for RCA tandem repeat detection; packaged as pip-installable Python package.

## Quick Start

```bash
# Install
conda activate ciderseq2
cd /path/to/ciderseq2
pip install .

# Run
ciderseq2 \
    --input sample.hifi.fasta \
    --genome reference_genome.fasta \
    --output results/ \
    --prefix my_run \
    --threads 8
```

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

### 3. Install ciderseq2

**Option A: Local install (recommended for servers without GitHub access)**

```bash
# On your local machine: clone the repo
git clone https://github.com/YxZhang-XHCY/ciderseq2.git
cd ciderseq2
git checkout tidehunter-deconcat

# Transfer to server
scp -r ciderseq2/ user@server:/path/to/tools/

# On the server: install
conda activate ciderseq2
cd /path/to/tools/ciderseq2
pip install .
```

**Option B: Direct install from GitHub**

```bash
pip install git+https://github.com/YxZhang-XHCY/ciderseq2.git@tidehunter-deconcat
```

### 4. Verify installation

```bash
ciderseq2 --help
ciderseq2 --version    # should print 3.0.0
```

## Usage

### Full pipeline (single sample)

```bash
ciderseq2 \
    --input sample.hifi.fasta \
    --genome reference_genome.fasta \
    --output results/ \
    --prefix my_run \
    --threads 8
```

### Multiple replicates

Comma-separated input files will be treated as biological replicates with automatic cross-replicate filtering:

```bash
ciderseq2 \
    --input rep1.fasta,rep2.fasta \
    --genome genome.fasta \
    --output results/ \
    --prefix multi_rep \
    --threads 8
```

### Skip DeConcat (use pre-computed results)

If you already have deconcat output files (`.fa` and `.stat`), use `--skip-deconcat`:

```bash
ciderseq2 \
    --input deconcat.fa:deconcat.stat \
    --genome genome.fasta \
    --output results/ \
    --prefix rerun \
    --skip-deconcat
```

### Standalone sub-commands

```bash
# Run TideHunter deconcat only
ciderseq2-deconcat \
    --input raw.hifi.fasta \
    --output-fa out.deconcat.fa \
    --output-stat out.deconcat.stat \
    --threads 8

# Generate BED file from existing results
ciderseq2-to-bed \
    --assess output/assess1.fasta \
    --blast output/reorder1.fasta_on_genomedb.blastn6 \
    --stat deconcat.stat \
    --output eccdna.bed
```

## Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `--input` | required | HiFi FASTA file(s), comma-separated for replicates |
| `--genome` | required | Host reference genome FASTA |
| `--output` | required | Output directory |
| `--prefix` | required | Output file prefix |
| `--threads` | 8 | Number of BLAST threads |
| `--gap-window` | 150 | Gap window for merging partial BLAST hits |
| `--tidehunter` | auto-detect | Path to TideHunter executable |
| `--skip-deconcat` | off | Skip DeConcat step; input format becomes `fasta:stat` |

## Output

```
<output>/<prefix>/
├── 01_deconcat/
│   ├── <sample>.deconcat.fa          # De-concatenated reads (consensus)
│   └── <sample>.deconcat.stat        # Copy number statistics (11-col TSV)
├── 02_eccDNA/output/
│   ├── <prefix>_eccDNA_candidates.fasta   ← Final eccDNA candidate sequences
│   ├── <prefix>_eccDNA.bed                ← Genomic coordinates (BED format)
│   ├── assess*.fasta                      # Per-sample assessed candidates
│   ├── reorder*.fasta                     # Reordered reads
│   └── rfile*.fasta                       # Filtered reads (copy > 0)
└── logs/
```

### Key output files

- **`<prefix>_eccDNA_candidates.fasta`** — Detected eccDNA consensus sequences
- **`<prefix>_eccDNA.bed`** — Genomic coordinates with 8 columns:

| Column | Description |
|--------|-------------|
| chrom | Chromosome |
| start | Start position (0-based) |
| end | End position |
| name | Read ID |
| score | Alignment identity × 10 |
| strand | + or - |
| length | eccDNA length (bp) |
| copies | RCA tandem repeat copy number |

## Dependencies

| Tool | Version | Install |
|------|---------|---------|
| Python | >=3.7 | conda |
| TideHunter | >=1.5 | `conda install -c bioconda tidehunter` |
| BLAST+ | >=2.10 | `conda install -c bioconda blast` |
| Biopython | <=1.77 | `pip install biopython==1.77` |
| Click | any | `pip install click` |

**Note:** Biopython must be version 1.77 or earlier (the `Bio.Blast.Applications` module was removed in later versions).

## Pipeline overview

1. **DeConcat (TideHunter)** — Detect RCA tandem repeats in HiFi reads; extract consensus sequences for circular reads
2. **eccDNA Detection (cs-eccDNA)** — BLAST against host genome, reorder partial hits, assess single-hit coverage, cross-replicate filtering
3. **BED Generation** — Map eccDNA candidates to genomic coordinates

## Citation

Mehta D, Cornet L, Hirsch-Hoffmann M, Zaidi, SSA, Vanderschuren H (2020) _Full-length sequencing of circular DNA viruses and extra-chromosomal circular DNA using CIDER-Seq._ **Nature Protocols**, Volume 15, pages 1673–1689; [https://doi.org/10.1038/s41596-020-0301-0](https://doi.org/10.1038/s41596-020-0301-0)

Mehta D, Hirsch-Hoffmann M, Were M, Patrignani A, Were H, Gruissem W, Vanderschuren H (2019) _A new full-length circular DNA sequencing method for viral-sized genomes reveals that RNAi transgenic plants provoke a shift in geminivirus populations in the field._ **Nucleic Acids Research**, Volume 47, Issue 2; [https://doi.org/10.1093/nar/gky914](https://doi.org/10.1093/nar/gky914)

## License

Copyright &copy; 2019, Matthias Hirsch-Hoffmann and Devang Mehta. Licensed under [GNU GPL v3](LICENSE.txt).
