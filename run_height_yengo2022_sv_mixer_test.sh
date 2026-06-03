#!/bin/bash
#SBATCH --job-name=sv_mixer_eur
#SBATCH --account=nn9114k
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err

set -euo pipefail

REPO_DIR="/cluster/projects/nn9114k/datngu/projects/ofrei_ldsc"
SUMSTATS="/cluster/projects/nn9114k/datngu/database/alkes_group_data/sumstats_107/PASS.Height.Yengo2022.sumstats.gz"
REF_DIR="${REPO_DIR}/ldsc_mixer_sv_reference"
OUT_DIR="${REPO_DIR}/results/height_yengo2022_sv_mixer_test"
OUT_PREFIX="${OUT_DIR}/PASS.Height.Yengo2022.sv_mixer_eur"

mkdir -p "${REPO_DIR}/logs" "${OUT_DIR}"
cd "${REPO_DIR}"

echo "Started: $(date)"
echo "Host: $(hostname)"
echo "Python:"
conda run -n ldsc python -c 'import sys; print(sys.executable); print(sys.version)'

conda run -n ldsc python "${REPO_DIR}/ldsc.py" \
  --h2 "${SUMSTATS}" \
  --ref-ld-chr "${REF_DIR}/ldscore/chr" \
  --w-ld-chr "${REF_DIR}/weights/chr" \
  --out "${OUT_PREFIX}"

echo "Finished: $(date)"
echo "LDSC log: ${OUT_PREFIX}.log"
