#!/bin/bash
#===============================================================================
# ciderseq2_pipeline.sh
# CIDER-Seq2 eccDNA Detection Pipeline
#
# Pipeline:
#   Step 1: DeConcat (TideHunter) - identify circular reads from RCA tandem repeats
#   Step 2: eccDNA detection (cs-eccDNA.py) - BLAST against host genome & filtering
#
# Usage:
#   bash ciderseq2_pipeline.sh \
#       -i sample1.fasta \
#       -r genome.fasta \
#       -o /path/to/output \
#       -p my_prefix \
#       -t 8
#
# Required:
#   -i  Input HiFi fasta files (comma-separated for multiple replicates)
#       Single sample:    -i sample.fasta
#       Multi replicate:  -i rep1.fasta,rep2.fasta
#   -r  Reference genome fasta
#   -o  Output directory
#   -p  Output prefix
#
# Optional:
#   -t  BLAST threads (default: 8)
#   -g  Gap window for merging partial BLAST hits (default: 150)
#   -s  Skip DeConcat step, input files are already deconcat outputs
#       In this mode, -i should be: fasta1:stat1[,fasta2:stat2,...]
#   -S  Path to CIDER-Seq2 installation (default: auto-detect from script location)
#   -T  Path to TideHunter executable (default: auto-detect from PATH)
#
# Output structure:
#   <outdir>/<prefix>/
#   ├── 01_deconcat/
#   │   ├── sample.deconcat.fa       # de-concatenated reads
#   │   └── sample.deconcat.stat     # rounds statistics (col2: rounds of concatenation)
#   ├── 02_eccDNA/
#   │   ├── output/                  # cs-eccDNA.py raw output
#   │   │   ├── rfile*.fasta         # filtered reads (rounds > 0)
#   │   │   ├── reorder*.fasta       # reordered reads
#   │   │   ├── assess*.fasta        # eccDNA candidates (final output for single sample)
#   │   │   └── replicatefile*.fasta # cross-replicate eccDNA (multi-sample only)
#   │   └── samples.list             # input list for cs-eccDNA.py
#   └── logs/
#       ├── deconcat_*.log
#       └── pipeline.log
#
# Environment: conda activate ciderseq2
#   (Python 3.7+, Biopython<=1.77, BLAST+, TideHunter, Click)
#===============================================================================

set -euo pipefail

#--- Default parameters --------------------------------------------------------
THREADS=8
GAP_WINDOW=150
SKIP_DECONCAT=false
CIDERSEQ_DIR=""
TIDEHUNTER_EXE=""

#--- Color output --------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info()  { echo -e "${GREEN}[INFO $(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }
log_warn()  { echo -e "${YELLOW}[WARN $(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }
log_error() { echo -e "${RED}[ERROR $(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }
log_step()  { echo -e "${BLUE}[STEP $(date '+%Y-%m-%d %H:%M:%S')]${NC} $*"; }

#--- Usage ---------------------------------------------------------------------
usage() {
    sed -n '/^# Usage:/,/^# Environment:/p' "$0" | sed 's/^# //' | sed 's/^#//'
    exit 1
}

#--- Parse arguments -----------------------------------------------------------
while getopts "i:r:o:p:t:g:sS:T:h" opt; do
    case $opt in
        i) INPUT_FILES="$OPTARG" ;;
        r) GENOME="$OPTARG" ;;
        o) OUTDIR="$OPTARG" ;;
        p) PREFIX="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        g) GAP_WINDOW="$OPTARG" ;;
        s) SKIP_DECONCAT=true ;;
        S) CIDERSEQ_DIR="$OPTARG" ;;
        T) TIDEHUNTER_EXE="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

#--- Validate required parameters ----------------------------------------------
if [[ -z "${INPUT_FILES:-}" || -z "${GENOME:-}" || -z "${OUTDIR:-}" || -z "${PREFIX:-}" ]]; then
    log_error "Missing required parameters: -i, -r, -o, -p are all required"
    usage
fi

#--- Auto-detect CIDER-Seq2 installation path ----------------------------------
if [[ -z "$CIDERSEQ_DIR" ]]; then
    # Try common locations
    SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
    for candidate in \
        "$SCRIPT_DIR" \
        "$SCRIPT_DIR/ciderseq2" \
        "$(dirname "$SCRIPT_DIR")" \
        "$(dirname "$SCRIPT_DIR")/ciderseq2" \
        "$HOME/tools/ciderseq2" \
        "$HOME/software/ciderseq2"; do
        if [[ -f "$candidate/tidehunter_deconcat.py" && -f "$candidate/eccDNA/cs-eccDNA.py" ]]; then
            CIDERSEQ_DIR="$candidate"
            break
        fi
    done
fi

if [[ -z "$CIDERSEQ_DIR" || ! -f "$CIDERSEQ_DIR/tidehunter_deconcat.py" ]]; then
    log_error "Cannot find CIDER-Seq2 installation. Use -S to specify path."
    exit 1
fi

CIDERSEQ_DIR="$(cd "$CIDERSEQ_DIR" && pwd)"
log_info "CIDER-Seq2 installation: $CIDERSEQ_DIR"

#--- Resolve absolute paths ----------------------------------------------------
GENOME="$(readlink -f "$GENOME")"
if [[ ! -f "$GENOME" ]]; then
    log_error "Genome file not found: $GENOME"
    exit 1
fi

#--- Setup directory structure -------------------------------------------------
WORKDIR="${OUTDIR}/${PREFIX}"
DIR_DECONCAT="${WORKDIR}/01_deconcat"
DIR_ECCDNA="${WORKDIR}/02_eccDNA"
DIR_LOG="${WORKDIR}/logs"
LOGFILE="${DIR_LOG}/pipeline.log"

mkdir -p "$DIR_DECONCAT" "$DIR_ECCDNA/output" "$DIR_LOG"

# Tee all output to log file
exec > >(tee -a "$LOGFILE") 2>&1

log_info "============================================================"
log_info "CIDER-Seq2 eccDNA Detection Pipeline"
log_info "============================================================"
log_info "Input files:    $INPUT_FILES"
log_info "Genome:         $GENOME"
log_info "Output:         $WORKDIR"
log_info "Prefix:         $PREFIX"
log_info "Threads:        $THREADS"
log_info "Gap window:     $GAP_WINDOW"
log_info "Skip DeConcat:  $SKIP_DECONCAT"
log_info "============================================================"

#--- Check dependencies --------------------------------------------------------
log_step "Checking dependencies..."
for cmd in python3 blastn makeblastdb; do
    if ! command -v "$cmd" &>/dev/null; then
        log_error "Required command not found: $cmd"
        log_error "Please activate conda environment: conda activate ciderseq2"
        exit 1
    fi
done

# Auto-detect TideHunter if not specified
if [[ -z "$TIDEHUNTER_EXE" ]]; then
    TIDEHUNTER_EXE="$(which TideHunter 2>/dev/null || true)"
    if [[ -z "$TIDEHUNTER_EXE" ]]; then
        log_error "TideHunter not found in PATH. Use -T to specify path."
        exit 1
    fi
fi
if [[ ! -x "$TIDEHUNTER_EXE" ]]; then
    log_error "TideHunter not executable: $TIDEHUNTER_EXE"
    exit 1
fi
log_info "TideHunter: $TIDEHUNTER_EXE"

# Check Biopython version
BIOPYTHON_VER=$(python3 -c "import Bio; print(Bio.__version__)" 2>/dev/null || echo "not_found")
if [[ "$BIOPYTHON_VER" == "not_found" ]]; then
    log_error "Biopython not installed"
    exit 1
fi
log_info "Biopython version: $BIOPYTHON_VER"

# Check Bio.Blast.Applications availability (removed in Biopython > 1.77)
if ! python3 -c "from Bio.Blast.Applications import NcbiblastnCommandline" 2>/dev/null; then
    log_error "Bio.Blast.Applications not available. Biopython must be <= 1.77"
    log_error "Current version: $BIOPYTHON_VER"
    exit 1
fi
log_info "Dependencies check passed"

#--- Parse input files ---------------------------------------------------------
IFS=',' read -ra SAMPLES <<< "$INPUT_FILES"
NUM_SAMPLES=${#SAMPLES[@]}
log_info "Number of samples: $NUM_SAMPLES"

if [[ $NUM_SAMPLES -eq 1 ]]; then
    log_info "Single sample mode: eccDNA candidates from assess file (no cross-replicate filtering)"
else
    log_info "Multi-sample mode: cross-replicate filtering enabled"
fi

#===============================================================================
# STEP 1: DeConcat - Identify circular reads (TideHunter)
#===============================================================================
# Prepare arrays for Step 2 input
declare -a DECONCAT_FASTAS
declare -a DECONCAT_STATS

if [[ "$SKIP_DECONCAT" == true ]]; then
    log_step "Skipping DeConcat (using pre-computed outputs)"
    # Parse format: fasta1:stat1,fasta2:stat2
    for sample in "${SAMPLES[@]}"; do
        IFS=':' read -r fasta stat <<< "$sample"
        fasta="$(readlink -f "$fasta")"
        stat="$(readlink -f "$stat")"
        if [[ ! -f "$fasta" ]]; then
            log_error "Fasta file not found: $fasta"
            exit 1
        fi
        if [[ ! -f "$stat" ]]; then
            log_error "Stat file not found: $stat"
            exit 1
        fi
        DECONCAT_FASTAS+=("$fasta")
        DECONCAT_STATS+=("$stat")
        log_info "  Sample: $fasta | $stat"
    done
else
    log_step "Step 1: DeConcat (TideHunter) - Identifying circular reads"

    SAMPLE_IDX=0
    for sample in "${SAMPLES[@]}"; do
        SAMPLE_IDX=$((SAMPLE_IDX + 1))
        SAMPLE_FILE="$(readlink -f "$sample")"
        SAMPLE_NAME="$(basename "${SAMPLE_FILE%.*}")"

        if [[ ! -f "$SAMPLE_FILE" ]]; then
            log_error "Input file not found: $SAMPLE_FILE"
            exit 1
        fi

        log_info "Processing sample ${SAMPLE_IDX}/${NUM_SAMPLES}: $SAMPLE_NAME"

        # Count input reads
        NUM_READS=$(grep -c "^>" "$SAMPLE_FILE" || true)
        log_info "  Input reads: $NUM_READS"

        # Output paths
        CLEAN_FA="${DIR_DECONCAT}/${SAMPLE_NAME}.deconcat.fa"
        CLEAN_STAT="${DIR_DECONCAT}/${SAMPLE_NAME}.deconcat.stat"

        # Run TideHunter-based DeConcat
        log_info "  Running TideHunter DeConcat..."

        START_TIME=$(date +%s)

        python3 "${CIDERSEQ_DIR}/tidehunter_deconcat.py" \
            --input "$SAMPLE_FILE" \
            --output-fa "$CLEAN_FA" \
            --output-stat "$CLEAN_STAT" \
            --tidehunter "$TIDEHUNTER_EXE" \
            --threads "$THREADS" \
            --min-copy 2 \
            --min-period 30 \
            2>&1 | tee "${DIR_LOG}/deconcat_${SAMPLE_NAME}.log"

        END_TIME=$(date +%s)
        ELAPSED=$((END_TIME - START_TIME))
        log_info "  DeConcat finished in ${ELAPSED}s"

        if [[ ! -f "$CLEAN_FA" ]]; then
            log_error "  DeConcat fasta not found: $CLEAN_FA"
            exit 1
        fi

        # Count output
        NUM_DECONCAT=$(grep -c "^>" "$CLEAN_FA" || true)
        NUM_CIRCULAR=$(awk -F'\t' '$2 > 0' "$CLEAN_STAT" | wc -l || true)
        log_info "  DeConcat reads: $NUM_DECONCAT (circular: $NUM_CIRCULAR)"

        DECONCAT_FASTAS+=("$CLEAN_FA")
        DECONCAT_STATS+=("$CLEAN_STAT")
    done
fi

log_info "DeConcat outputs prepared: ${#DECONCAT_FASTAS[@]} samples"

#===============================================================================
# STEP 2: eccDNA Detection - BLAST against host genome & cross-replicate analysis
#===============================================================================
log_step "Step 2: eccDNA Detection"

cd "$DIR_ECCDNA"

# Create list file (tab-separated: fasta<TAB>stat)
LIST_FILE="${DIR_ECCDNA}/samples.list"
> "$LIST_FILE"
for i in "${!DECONCAT_FASTAS[@]}"; do
    echo -e "${DECONCAT_FASTAS[$i]}\t${DECONCAT_STATS[$i]}" >> "$LIST_FILE"
done

log_info "List file created: $LIST_FILE"
cat "$LIST_FILE"

# Run cs-eccDNA.py
log_info "Running cs-eccDNA.py (threads=${THREADS}, gap_window=${GAP_WINDOW})..."
START_TIME=$(date +%s)

python3 "${CIDERSEQ_DIR}/eccDNA/cs-eccDNA.py" \
    "$LIST_FILE" \
    "$GENOME" \
    --blast_threads "$THREADS" \
    --gap_window "$GAP_WINDOW" \
    2>&1 | tee "${DIR_LOG}/cs-eccDNA.log"

END_TIME=$(date +%s)
ELAPSED=$((END_TIME - START_TIME))
log_info "cs-eccDNA.py finished in ${ELAPSED}s"

#===============================================================================
# STEP 3: Summarize results
#===============================================================================
log_step "Summarizing results"

cd "$DIR_ECCDNA"

echo ""
log_info "============================================================"
log_info "Pipeline Summary"
log_info "============================================================"

# Per-sample intermediate results
for i in $(seq 1 ${#DECONCAT_FASTAS[@]}); do
    RFILE="output/rfile${i}.fasta"
    ASSESS="output/assess${i}.fasta"
    if [[ -f "$RFILE" ]]; then
        N_RFILE=$(grep -c "^>" "$RFILE" 2>/dev/null || echo 0)
        log_info "  Sample ${i} - filtered reads (rounds>0): $N_RFILE"
    fi
    if [[ -f "$ASSESS" ]]; then
        N_ASSESS=$(grep -c "^>" "$ASSESS" 2>/dev/null || echo 0)
        log_info "  Sample ${i} - assessed eccDNA candidates: $N_ASSESS"
    fi
done

# Final replicate files (multi-sample) or assess files (single-sample)
if [[ $NUM_SAMPLES -ge 2 ]]; then
    FINAL_FILES=($(ls output/replicatefile*.fasta 2>/dev/null || true))
    if [[ ${#FINAL_FILES[@]} -gt 0 ]]; then
        for f in "${FINAL_FILES[@]}"; do
            N=$(grep -c "^>" "$f" 2>/dev/null || echo 0)
            log_info "  Cross-replicate eccDNA ($(basename "$f")): $N"
        done
    else
        log_warn "Cross-replicate comparison produced no shared eccDNA"
        log_info "  Per-sample candidates in: output/assess*.fasta"
    fi
else
    # Single sample: assess1.fasta is the final result
    ASSESS="output/assess1.fasta"
    if [[ -f "$ASSESS" ]]; then
        N_ASSESS=$(grep -c "^>" "$ASSESS" 2>/dev/null || echo 0)
        log_info "  Final eccDNA candidates (single sample): $N_ASSESS"
        # Copy to a cleaner final output name
        FINAL_OUT="output/${PREFIX}_eccDNA_candidates.fasta"
        cp "$ASSESS" "$FINAL_OUT"
        log_info "  Saved as: $FINAL_OUT"
    else
        log_warn "No eccDNA candidates detected"
    fi
fi

#--- Generate BED file with genomic coordinates --------------------------------
log_step "Generating eccDNA BED file"

BED_OUT="output/${PREFIX}_eccDNA.bed"
BED_ASSESS_ARGS=""
BED_BLAST_ARGS=""
BED_STAT_ARGS=""

for i in $(seq 1 ${#DECONCAT_FASTAS[@]}); do
    ASSESS="output/assess${i}.fasta"
    BLAST="output/reorder${i}.fasta_on_genomedb.blastn6"
    if [[ -f "$ASSESS" && -f "$BLAST" ]]; then
        BED_ASSESS_ARGS+=" --assess $ASSESS"
        BED_BLAST_ARGS+=" --blast $BLAST"
    fi
    if [[ -f "${DECONCAT_STATS[$i-1]}" ]]; then
        BED_STAT_ARGS+=" --stat ${DECONCAT_STATS[$i-1]}"
    fi
done

if [[ -n "$BED_ASSESS_ARGS" ]]; then
    python3 "${CIDERSEQ_DIR}/eccDNA_to_bed.py" \
        $BED_ASSESS_ARGS \
        $BED_BLAST_ARGS \
        $BED_STAT_ARGS \
        --output "$BED_OUT" 2>&1
    if [[ -f "$BED_OUT" ]]; then
        N_BED=$(wc -l < "$BED_OUT" | tr -d ' ')
        log_info "  BED file: $BED_OUT ($N_BED eccDNA regions)"
    fi
fi

echo ""
log_info "============================================================"
log_info "Results directory: $DIR_ECCDNA/output/"
log_info "Logs directory:    $DIR_LOG/"
log_info "============================================================"

# List output files
log_info "Output files:"
ls -lh "$DIR_ECCDNA/output/"*.fasta "$DIR_ECCDNA/output/"*.bed 2>/dev/null | while read line; do
    log_info "  $line"
done

log_info "Pipeline complete!"
