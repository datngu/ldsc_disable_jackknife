#!/bin/bash

set -euo pipefail

LD_DIR="/Users/thanhdng/github/mixer_sv/1kg_combined_plink_eur_AC3"
ANNOT="${LD_DIR}/annot_mat.txt"
SNP_FILE="${LD_DIR}/SV_variants.txt"
OUT_DIR="ldsc_mixer_sv_reference"

python make_ldsc_annot.py \
  --ld-dir "${LD_DIR}" \
  --annot "${ANNOT}" \
  --snp-file "${SNP_FILE}" \
  --roi-mode replace \
  --out-dir "${OUT_DIR}" \
  --chromosomes 1-22

echo "LDSC reference files written to: ${OUT_DIR}"
echo "Use with LDSC:"
echo "  --ref-ld-chr ${OUT_DIR}/ldscore/chr"
echo "  --w-ld-chr   ${OUT_DIR}/weights/chr"
