#!/usr/bin/env python3
"""
Create LDSC-compatible annotation, LD-score, and weight files from MiXER-SV LD.

LDSC annotation LD scores are computed with the same LD matrix and annotation
matrix used by MiXER-SV:
    LD_score_k = R^2 @ annotation_k

The weight file uses the same LD score that MiXER-SV uses for inverse-LD
weights when all SNPs are defined:
    mixer_weight = 1 / (R^2 @ 1)

By default, chr6:25,000,000-34,000,000 is excluded to mirror MiXER-SV's
MHC filter.

Outputs are written per chromosome so they can be passed to LDSC with
--ref-ld-chr and --w-ld-chr prefixes.
"""

import argparse
import gzip
import sys
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

sys.path.insert(0, ".")
from utils.ld_matrix import LDMatrix


def parse_chromosomes(value: str) -> list[int]:
    chroms: list[int] = []
    for part in value.split(","):
        part = part.strip()
        if not part:
            continue
        if "-" in part:
            start, end = part.split("-", 1)
            chroms.extend(range(int(start), int(end) + 1))
        else:
            chroms.append(int(part))
    return sorted(set(chroms))


def load_annotations(path: Path, only_base: bool, snp_file: Path | None,
                     roi_mode: str) -> pd.DataFrame:
    annot = pd.read_csv(path, sep="\t")
    if "snp" not in annot.columns:
        raise ValueError(f"{path} must contain a 'snp' column")
    annot = annot.set_index("snp")
    if only_base:
        annot = annot.iloc[:, [0]]
    if snp_file is not None:
        roi_snps = set(pd.read_csv(snp_file, header=None)[0].astype(str).values)
        within = annot.index.astype(str).isin(roi_snps).astype(int)
        roi_annot = pd.DataFrame({
            "outside_region": 1 - within,
            "within_region": within,
        }, index=annot.index)
        if roi_mode == "replace":
            annot = roi_annot
        elif roi_mode == "append":
            annot = pd.concat([annot, roi_annot], axis=1)
        else:
            raise ValueError(f"Unknown roi_mode: {roi_mode}")
    return annot


def align_annotations(info: pd.DataFrame, annot: pd.DataFrame) -> pd.DataFrame:
    missing = ~info["snp"].isin(annot.index)
    if missing.any():
        examples = ", ".join(info.loc[missing, "snp"].head(5).astype(str))
        raise ValueError(
            f"{missing.sum():,} LD SNPs are missing from the annotation matrix. "
            f"Examples: {examples}"
        )
    return annot.loc[info["snp"]]


def write_gzip_tsv(df: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", newline="") as f:
        df.to_csv(f, sep="\t", index=False, float_format="%.10g")


def write_vector(values: Iterable[float], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.asarray(list(values), dtype=float)
    np.savetxt(path, arr.reshape(1, -1), fmt="%.10g", delimiter="\t")


def keep_mask_after_region_filter(
    info: pd.DataFrame,
    exclude_mhc: bool,
    mhc_chr: str,
    mhc_start: int,
    mhc_end: int,
) -> np.ndarray:
    if not exclude_mhc:
        return np.ones(len(info), dtype=bool)

    chrom = info["chr"].astype(str).str.removeprefix("chr").str.removeprefix("CHR")
    in_mhc = (
        (chrom == str(mhc_chr).removeprefix("chr").removeprefix("CHR")) &
        info["pos"].between(mhc_start, mhc_end)
    )
    return ~in_mhc.values


def write_chromosome_outputs(
    ld_dir: Path,
    annot: pd.DataFrame,
    chrom: int,
    output_dir: Path,
    exclude_mhc: bool,
    mhc_chr: str,
    mhc_start: int,
    mhc_end: int,
) -> dict[str, object]:
    chr_ld_dir = ld_dir / f"chr{chrom}.ldmat"
    if not chr_ld_dir.exists():
        raise FileNotFoundError(f"Missing LD matrix directory: {chr_ld_dir}")

    print(f"Loading chr{chrom}: {chr_ld_dir}")
    ld_obj = LDMatrix.load_data(chr_ld_dir)
    info = ld_obj.info.copy()
    chr_annot = align_annotations(info, annot)
    keep_mask = keep_mask_after_region_filter(
        info=info,
        exclude_mhc=exclude_mhc,
        mhc_chr=mhc_chr,
        mhc_start=mhc_start,
        mhc_end=mhc_end,
    )

    annot_cols = list(chr_annot.columns)
    annot_values = chr_annot.values.astype(float)
    annot_values_for_scores = annot_values.copy()
    annot_values_for_scores[~keep_mask, :] = 0.0

    ld_r_sq = ld_obj.get_full_ld_matrix_sq()
    l2_values = np.asarray(ld_r_sq.dot(annot_values_for_scores))
    l2_cols = [f"{col}L2" for col in annot_cols]
    weight_l2 = np.asarray(ld_r_sq.dot(keep_mask.astype(int))).ravel()
    inv_weight = np.zeros_like(weight_l2, dtype=float)
    valid_weight = weight_l2 > 0
    inv_weight[valid_weight] = 1.0 / weight_l2[valid_weight]

    info_out = info.loc[keep_mask].reset_index(drop=True)
    chr_annot_out = chr_annot.loc[keep_mask].reset_index(drop=True)
    l2_values_out = l2_values[keep_mask, :]
    weight_l2_out = weight_l2[keep_mask]
    inv_weight_out = inv_weight[keep_mask]

    annot_df = pd.DataFrame({
        "CHR": info_out["chr"].astype(str).values,
        "BP": info_out["pos"].values,
        "SNP": info_out["snp"].values,
        "CM": np.zeros(len(info_out), dtype=float),
    })
    for col in annot_cols:
        annot_df[col] = chr_annot_out[col].values

    ldscore_df = pd.DataFrame({
        "CHR": info_out["chr"].astype(str).values,
        "SNP": info_out["snp"].values,
        "BP": info_out["pos"].values,
    })
    for col, values in zip(l2_cols, l2_values_out.T):
        ldscore_df[col] = values

    weight_df = pd.DataFrame({
        "CHR": info_out["chr"].astype(str).values,
        "SNP": info_out["snp"].values,
        "BP": info_out["pos"].values,
        "L2": weight_l2_out,
    })

    mixer_weight_df = weight_df.copy()
    mixer_weight_df["INV_LD_WEIGHT"] = inv_weight_out

    write_gzip_tsv(annot_df, output_dir / "annot" / f"chr{chrom}.annot.gz")
    write_gzip_tsv(ldscore_df, output_dir / "ldscore" / f"chr{chrom}.l2.ldscore.gz")
    write_gzip_tsv(weight_df, output_dir / "weights" / f"chr{chrom}.l2.ldscore.gz")
    write_gzip_tsv(
        mixer_weight_df,
        output_dir / "weights" / f"chr{chrom}.mixer_inverse_ld_weights.tsv.gz",
    )

    annot_values_out = chr_annot_out.values.astype(float)
    m_values = annot_values_out.sum(axis=0)
    m_5_50_values = annot_values_out[info_out["maf"].values >= 0.05].sum(axis=0)

    write_vector(m_values, output_dir / "ldscore" / f"chr{chrom}.l2.M")
    write_vector(m_5_50_values, output_dir / "ldscore" / f"chr{chrom}.l2.M_5_50")

    raw_counts = annot_values_out.sum(axis=0)
    raw_counts_5_50 = annot_values_out[info_out["maf"].values >= 0.05].sum(axis=0)
    summary = pd.DataFrame({
        "annotation": annot_cols,
        "M": m_values,
        "M_5_50": m_5_50_values,
        "raw_annotation_sum": raw_counts,
        "raw_annotation_sum_5_50": raw_counts_5_50,
    })
    summary.to_csv(output_dir / "ldscore" / f"chr{chrom}.annotation_counts.tsv",
                   sep="\t", index=False, float_format="%.10g")

    return {
        "chrom": chrom,
        "n_snps": len(info_out),
        "n_snps_excluded": int((~keep_mask).sum()),
        "annotations": annot_cols,
        "mean_weight_ldscore": float(np.mean(weight_l2_out)),
        "min_weight_ldscore": float(np.min(weight_l2_out)),
        "max_weight_ldscore": float(np.max(weight_l2_out)),
    }


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Create LDSC files from MiXER-SV LD matrices and annotations."
    )
    p.add_argument("--ld-dir", type=Path, required=True,
                   help="Directory containing chr*.ldmat directories")
    p.add_argument("--annot", type=Path, required=True,
                   help="MiXER-SV annot_mat.txt with a snp column")
    p.add_argument("--snp-file", type=Path, default=None,
                   help="Optional MiXER-SV SV_variants.txt; adds within/outside ROI annotations")
    p.add_argument("--roi-mode", choices=["replace", "append"], default="replace",
                   help="With --snp-file, replace annotations with outside/within ROI columns or append them")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--chromosomes", default="1-22",
                   help="Chromosomes, e.g. '10' or '1-22' or '1,3,10'")
    p.add_argument("--only-base", action="store_true",
                   help="Use only the first annotation column")
    p.add_argument("--include-mhc", action="store_true",
                   help="Do not exclude chr6 MHC rows from LDSC outputs")
    p.add_argument("--mhc-chr", default="6")
    p.add_argument("--mhc-start", type=int, default=25_000_000)
    p.add_argument("--mhc-end", type=int, default=34_000_000)
    return p.parse_args()


def write_summary(summaries: list[dict[str, object]], output_dir: Path) -> None:
    summary_path = output_dir / "make_ldsc_annot_summary.tsv"
    new_summary = pd.DataFrame(summaries)
    if summary_path.exists():
        old_summary = pd.read_csv(summary_path, sep="\t")
        old_summary = old_summary.loc[~old_summary["chrom"].isin(new_summary["chrom"])]
        new_summary = pd.concat([old_summary, new_summary], ignore_index=True)
    new_summary = new_summary.sort_values("chrom")
    new_summary.to_csv(summary_path, sep="\t", index=False)


def main() -> None:
    args = parse_args()
    chroms = parse_chromosomes(args.chromosomes)

    print(f"Loading annotations: {args.annot}")
    annot = load_annotations(args.annot, args.only_base, args.snp_file, args.roi_mode)
    print(f"Annotation columns: {list(annot.columns)}")
    print("LD-score mode: plain LDSC R^2 @ annotation")
    if args.include_mhc:
        print("MHC filter: disabled")
    else:
        print(f"MHC filter: excluding chr{args.mhc_chr}:{args.mhc_start}-{args.mhc_end}")

    summaries = []
    for chrom in chroms:
        summaries.append(
            write_chromosome_outputs(
                ld_dir=args.ld_dir,
                annot=annot,
                chrom=chrom,
                output_dir=args.out_dir,
                exclude_mhc=not args.include_mhc,
                mhc_chr=args.mhc_chr,
                mhc_start=args.mhc_start,
                mhc_end=args.mhc_end,
            )
        )

    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_summary(summaries, args.out_dir)

    print(f"\nDone. Outputs written to: {args.out_dir}")
    print("Use these LDSC prefixes:")
    print(f"  --ref-ld-chr {args.out_dir / 'ldscore' / 'chr'}")
    print(f"  --w-ld-chr   {args.out_dir / 'weights' / 'chr'}")


if __name__ == "__main__":
    main()
