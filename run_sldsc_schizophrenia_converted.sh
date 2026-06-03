#!/bin/bash
#SBATCH --job-name=sldsc_scz_sv
#SBATCH --account=nn9114k
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO_DIR="/cluster/projects/nn9114k/datngu/projects/ofrei_ldsc"
SUMSTATS="${REPO_DIR}/results/sldsc_sumstats_107/FIXED.Schizophrenia.Trubetskoy2022.converted.sumstats.gz"
REF_DIR="${REPO_DIR}/ldsc_mixer_sv_reference"
OUT_DIR="${REPO_DIR}/results/sldsc_sumstats_107"
OUT_PREFIX="${OUT_DIR}/FIXED.Schizophrenia.Trubetskoy2022.sldsc_sv"

mkdir -p "${REPO_DIR}/logs" "${OUT_DIR}"
cd "${REPO_DIR}"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Converted sumstats: ${SUMSTATS}"

conda run -n ldsc python "${REPO_DIR}/ldsc.py" \
  --h2 "${SUMSTATS}" \
  --ref-ld-chr "${REF_DIR}/ldscore/chr" \
  --w-ld-chr "${REF_DIR}/weights/chr" \
  --out "${OUT_PREFIX}"

echo "Finished: $(date)"
echo "LDSC log: ${OUT_PREFIX}.log"
