#!/bin/bash -l
#$ -N fly_null_bins
#$ -cwd
#$ -j y
#$ -pe omp 4
#$ -l h_rt=12:00:00

set -euo pipefail

MINICONDA_MODULE="${MINICONDA_MODULE:-miniconda}"
CONDA_ENV_NAME="${CONDA_ENV_NAME:-shadow_null_env}"
PYTHON_BIN="${PYTHON_BIN:-python}"

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
FLY_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
PAPER_ROOT="$(cd "$FLY_DIR/../../.." && pwd)"
PYTHON_SCRIPT="${PYTHON_SCRIPT:-$FLY_DIR/original_scripts/make_null_enhancer_props_v2_fly_direct_from_flilesoutput.py}"

ENHANCER_ROOT="${ENHANCER_ROOT:-$FLY_DIR/data/enhancer_pairs}"
FLANK_MANIFEST="${FLANK_MANIFEST:-$PAPER_ROOT/duplicaton_analysis/FLY_flank_analysis/shadow_pairs_flanks2/manifest.tsv}"
FLANKS_ROOT="${FLANKS_ROOT:-$PAPER_ROOT/duplicaton_analysis/FLY_flank_analysis/shadow_pairs_flanks2}"
OUT_DIR="${OUT_DIR:-$FLY_DIR/outputs/null_bin_runs}"

OBSERVED_PAIRS_CSV="${OBSERVED_PAIRS_CSV:-$FLY_DIR/data/FINAL_breakdown_single_double_flanks_UPDATED_9.9.e-5_flankmethod.csv}"
GENE_SIZES_CSV="${GENE_SIZES_CSV:-}"

REPS="${REPS:-100}"
REP_WORKERS="${REP_WORKERS:-${NSLOTS:-1}}"
BUCKETS="${BUCKETS:-2,3,>=4}"
SEED="${SEED:-42}"
ENHANCER_HIT_EVALUE="${ENHANCER_HIT_EVALUE:-9.9e-5}"
DOUBLE_FLANK_HIT_EVALUE="${DOUBLE_FLANK_HIT_EVALUE:-9.9e-5}"

SKIP_OBSERVED="${SKIP_OBSERVED:-0}"
NO_PLOT="${NO_PLOT:-1}"
VERBOSE="${VERBOSE:-1}"
DRY_RUN="${DRY_RUN:-0}"

module purge
module load "$MINICONDA_MODULE"

if ! command -v conda >/dev/null 2>&1; then
  echo "[ERROR] conda is not available after loading module: $MINICONDA_MODULE" >&2
  exit 1
fi

CONDA_BASE="$(conda info --base)"
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV_NAME"

echo "Job started on $(hostname) at $(date)"
echo "Python script: $PYTHON_SCRIPT"
echo "Enhancer root: $ENHANCER_ROOT"
echo "Flank manifest: $FLANK_MANIFEST"
echo "Flanks root: $FLANKS_ROOT"
echo "Observed pairs CSV: $OBSERVED_PAIRS_CSV"
echo "Out dir: $OUT_DIR"
echo "Enhancer hit evalue: $ENHANCER_HIT_EVALUE"
echo "Double flank hit evalue: $DOUBLE_FLANK_HIT_EVALUE"
echo "Buckets: $BUCKETS"
echo "Reps: $REPS"
echo "Rep workers: $REP_WORKERS"

if [[ -z "$ENHANCER_ROOT" || -z "$FLANK_MANIFEST" || -z "$FLANKS_ROOT" || -z "$OUT_DIR" ]]; then
  echo "[ERROR] Required variables: ENHANCER_ROOT, FLANK_MANIFEST, FLANKS_ROOT, OUT_DIR" >&2
  exit 1
fi

if [[ ! -f "$PYTHON_SCRIPT" ]]; then
  echo "[ERROR] Python script not found: $PYTHON_SCRIPT" >&2
  exit 1
fi

if [[ ! -d "$ENHANCER_ROOT" ]]; then
  echo "[ERROR] Enhancer root not found: $ENHANCER_ROOT" >&2
  exit 1
fi

if [[ ! -f "$FLANK_MANIFEST" ]]; then
  echo "[ERROR] Flank manifest not found: $FLANK_MANIFEST" >&2
  exit 1
fi

if [[ ! -d "$FLANKS_ROOT" ]]; then
  echo "[ERROR] Flanks root not found: $FLANKS_ROOT" >&2
  exit 1
fi

if [[ "$SKIP_OBSERVED" == "0" && ! -f "$OBSERVED_PAIRS_CSV" ]]; then
  echo "[ERROR] Observed pairs CSV not found: $OBSERVED_PAIRS_CSV" >&2
  exit 1
fi

if ! command -v "$PYTHON_BIN" >/dev/null 2>&1; then
  echo "[ERROR] Python executable not found on PATH: $PYTHON_BIN" >&2
  exit 1
fi

if ! command -v blastn >/dev/null 2>&1; then
  echo "[ERROR] blastn is not on PATH in env: $CONDA_ENV_NAME" >&2
  exit 1
fi

mkdir -p "$OUT_DIR"

cmd=(
  "$PYTHON_BIN" "$PYTHON_SCRIPT"
  --enhancer_root "$ENHANCER_ROOT"
  --flank_manifest "$FLANK_MANIFEST"
  --flanks_root "$FLANKS_ROOT"
  --out_dir "$OUT_DIR"
  --buckets "$BUCKETS"
  --reps "$REPS"
  --rep_workers "$REP_WORKERS"
  --seed "$SEED"
  --enhancer_hit_evalue "$ENHANCER_HIT_EVALUE"
  --double_flank_hit_evalue "$DOUBLE_FLANK_HIT_EVALUE"
)

if [[ "$SKIP_OBSERVED" == "1" ]]; then
  cmd+=(--skip_observed)
else
  cmd+=(--observed_pairs_csv "$OBSERVED_PAIRS_CSV")
fi

if [[ -n "$GENE_SIZES_CSV" ]]; then
  cmd+=(--gene_sizes_csv "$GENE_SIZES_CSV")
fi

if [[ "$NO_PLOT" == "1" ]]; then
  cmd+=(--no_plot)
fi

if [[ "$VERBOSE" == "1" ]]; then
  cmd+=(--verbose)
fi

if [[ "$DRY_RUN" == "1" ]]; then
  cmd+=(--dry_run)
fi

printf 'Running command:\n'
printf ' %q' "${cmd[@]}"
printf '\n'

"${cmd[@]}"

echo "Job finished at $(date)"
