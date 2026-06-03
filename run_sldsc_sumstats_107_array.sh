#!/bin/bash
#SBATCH --job-name=sldsc_sv107
#SBATCH --account=nn9114k
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-106%30
#SBATCH --output=logs/%x_%A_%a.out
#SBATCH --error=logs/%x_%A_%a.err

set -euo pipefail

REPO_DIR="/cluster/projects/nn9114k/datngu/projects/ofrei_ldsc"
MANIFEST="${REPO_DIR}/manifests/sumstats_107.txt"
REF_DIR="${REPO_DIR}/ldsc_mixer_sv_reference"
OUT_DIR="${REPO_DIR}/results/sldsc_sumstats_107"

mkdir -p "${REPO_DIR}/logs" "${OUT_DIR}"
cd "${REPO_DIR}"

SUMSTATS="$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${MANIFEST}")"
if [[ -z "${SUMSTATS}" ]]; then
  echo "No sumstats path for SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}" >&2
  exit 2
fi

TRAIT="$(basename "${SUMSTATS}" .sumstats.gz)"
OUT_PREFIX="${OUT_DIR}/${TRAIT}.sldsc_sv"

echo "Started: $(date)"
echo "Job: ${SLURM_JOB_ID:-NA}; task: ${SLURM_ARRAY_TASK_ID}"
echo "Host: $(hostname)"
echo "Trait: ${TRAIT}"
echo "Sumstats: ${SUMSTATS}"

conda run -n ldsc python "${REPO_DIR}/ldsc.py" \
  --h2 "${SUMSTATS}" \
  --ref-ld-chr "${REF_DIR}/ldscore/chr" \
  --w-ld-chr "${REF_DIR}/weights/chr" \
  --out "${OUT_PREFIX}"

echo "Finished: $(date)"
echo "LDSC log: ${OUT_PREFIX}.log"
