#!/usr/bin/env python
import argparse
import math
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


def parse_mixer_list(path):
    values = {}
    with open(path) as handle:
        for line in handle:
            if ":" not in line:
                continue
            key, value = line.strip().split(":", 1)
            value = value.strip()
            try:
                values[key] = float(value)
            except ValueError:
                values[key] = value
    return values


def parse_sldsc_log(path):
    text = path.read_text()
    values = {}
    patterns = {
        "ldsc_total_h2": r"Total Observed scale h2:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
        "ldsc_lambda_gc": r"Lambda GC:\s+([-+0-9.eE]+)",
        "ldsc_mean_chi2": r"Mean Chi\^2:\s+([-+0-9.eE]+)",
        "ldsc_intercept": r"Intercept:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
        "ldsc_ratio": r"Ratio:\s+([-+0-9.eE]+)\s+\(([-+0-9.eE]+)\)",
    }
    for key, pattern in patterns.items():
        match = re.search(pattern, text)
        if match:
            values[key] = float(match.group(1))
            if len(match.groups()) > 1:
                values[key + "_se"] = float(match.group(2))

    for label, prefix in [
        ("Observed scale h2", "ldsc_h2"),
        ("Observed scale h2 SE", "ldsc_h2_se"),
        ("Proportion of SNPs", "ldsc_prop_snps"),
        ("Proportion of h2g", "ldsc_prop_h2"),
        ("Enrichment", "ldsc_enrichment"),
        ("Coefficients", "ldsc_coef"),
        ("Coefficient SE", "ldsc_coef_se"),
    ]:
        match = re.search(rf"{re.escape(label)}:\s+([-+0-9.eE]+)\s+([-+0-9.eE]+)", text)
        if match:
            values[prefix + "_outside_region"] = float(match.group(1))
            values[prefix + "_within_region"] = float(match.group(2))

    return values


def trait_from_sumstats(path):
    name = Path(path).name
    suffix = ".sumstats.gz"
    return name[:-len(suffix)] if name.endswith(suffix) else name


def pearson(x, y):
    pair = pd.concat([x, y], axis=1).dropna()
    if len(pair) < 2:
        return math.nan
    return pair.iloc[:, 0].corr(pair.iloc[:, 1], method="pearson")


def spearman(x, y):
    pair = pd.concat([x, y], axis=1).dropna()
    if len(pair) < 2:
        return math.nan
    return pair.iloc[:, 0].corr(pair.iloc[:, 1], method="spearman")


def scatter(df, x_col, y_col, out_png, title, xlabel, ylabel):
    plot_df = df[[x_col, y_col, "trait"]].dropna()
    fig, ax = plt.subplots(figsize=(6, 5), dpi=160)
    ax.scatter(plot_df[x_col], plot_df[y_col], s=22, alpha=0.75)
    if len(plot_df) > 1:
        lo = min(plot_df[x_col].min(), plot_df[y_col].min())
        hi = max(plot_df[x_col].max(), plot_df[y_col].max())
        ax.plot([lo, hi], [lo, hi], color="black", linewidth=1, alpha=0.5)
        ax.text(
            0.03,
            0.97,
            f"n={len(plot_df)}\nPearson r={pearson(plot_df[x_col], plot_df[y_col]):.3f}\n"
            f"Spearman rho={spearman(plot_df[x_col], plot_df[y_col]):.3f}",
            transform=ax.transAxes,
            va="top",
            ha="left",
            fontsize=9,
            bbox={"facecolor": "white", "edgecolor": "0.85", "alpha": 0.9},
        )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.tight_layout()
    fig.savefig(out_png)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True)
    parser.add_argument("--ldsc-dir", required=True)
    parser.add_argument("--mixer-dir", required=True)
    parser.add_argument("--out-dir", required=True)
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = []
    with open(args.manifest) as handle:
        for line in handle:
            sumstats = line.strip()
            if not sumstats:
                continue
            trait = trait_from_sumstats(sumstats)
            row = {"trait": trait, "sumstats": sumstats}
            ldsc_log = Path(args.ldsc_dir) / f"{trait}.sldsc_sv.log"
            mixer_list = Path(args.mixer_dir) / f"univar_{trait}_list.txt"
            row["ldsc_log"] = str(ldsc_log)
            row["mixer_list"] = str(mixer_list)
            row["ldsc_status"] = "ok" if ldsc_log.exists() else "missing"
            row["mixer_status"] = "ok" if mixer_list.exists() else "missing"
            if ldsc_log.exists():
                row.update(parse_sldsc_log(ldsc_log))
            if mixer_list.exists():
                mixer = parse_mixer_list(mixer_list)
                row["mixer_enrichment_outside_region"] = mixer.get("enrich_vs_base_outside_region_trait1")
                row["mixer_enrichment_outside_region_se"] = mixer.get("enrich_vs_base_outside_region_trait1_se")
                row["mixer_enrichment_within_region"] = mixer.get("enrich_vs_base_within_region_trait1")
                row["mixer_enrichment_within_region_se"] = mixer.get("enrich_vs_base_within_region_trait1_se")
                row["mixer_total_h2"] = mixer.get("total_h2_trait1")
                row["mixer_total_h2_se"] = mixer.get("total_h2_trait1_se")
                row["mixer_h2_outside_region"] = mixer.get("h2_outside_region_trait1")
                row["mixer_h2_within_region"] = mixer.get("h2_within_region_trait1")
            rows.append(row)

    df = pd.DataFrame(rows)
    df.to_csv(out_dir / "sldsc_mixer_sv_comparison.tsv", sep="\t", index=False)

    complete = df[(df["ldsc_status"] == "ok") & (df["mixer_status"] == "ok")].copy()
    metrics = []
    for ldsc_col, mixer_col, label in [
        ("ldsc_enrichment_within_region", "mixer_enrichment_within_region", "within_region_enrichment"),
        ("ldsc_enrichment_outside_region", "mixer_enrichment_outside_region", "outside_region_enrichment"),
        ("ldsc_total_h2", "mixer_total_h2", "total_h2"),
    ]:
        metrics.append({
            "metric": label,
            "n": int(complete[[ldsc_col, mixer_col]].dropna().shape[0]),
            "pearson_r": pearson(complete[ldsc_col], complete[mixer_col]),
            "spearman_rho": spearman(complete[ldsc_col], complete[mixer_col]),
        })
    pd.DataFrame(metrics).to_csv(out_dir / "sldsc_mixer_sv_correlations.tsv", sep="\t", index=False)

    scatter(
        complete,
        "mixer_enrichment_within_region",
        "ldsc_enrichment_within_region",
        out_dir / "scatter_within_region_enrichment.png",
        "SV-region enrichment: MiXeR-SV vs sLDSC",
        "MiXeR-SV enrichment",
        "sLDSC enrichment",
    )
    scatter(
        complete,
        "mixer_total_h2",
        "ldsc_total_h2",
        out_dir / "scatter_total_h2.png",
        "Total h2: MiXeR-SV vs sLDSC",
        "MiXeR-SV total h2",
        "sLDSC total h2",
    )

    print(f"Wrote {out_dir / 'sldsc_mixer_sv_comparison.tsv'}")
    print(f"Wrote {out_dir / 'sldsc_mixer_sv_correlations.tsv'}")


if __name__ == "__main__":
    main()
