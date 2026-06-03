# LDSC Python 3 Fork

This repository contains a Python 3 port of LDSC v1.0.1 with local workflows for
stratified LD Score regression (sLDSC) against structural variant (SV) annotations.

The code has been tested in the `ldsc` conda environment on Saga with:

- Python 3.8.20
- NumPy 1.24.3
- pandas 1.5.3
- SciPy 1.10.1
- bitarray 2.6.0

## Setup

Create or update the conda environment:

```bash
conda env create --file environment.yml
conda activate ldsc
```

For an existing environment:

```bash
conda env update --file environment.yml
conda activate ldsc
```

Basic checks:

```bash
python ldsc.py -h
python munge_sumstats.py -h
nosetests -v test
```

The current Python 3 port passes the full test suite:

```text
Ran 167 tests
OK
```

## Standard LDSC Usage

Estimate SNP heritability from LDSC-formatted summary statistics:

```bash
python ldsc.py \
  --h2 trait.sumstats.gz \
  --ref-ld-chr path/to/ref/chr \
  --w-ld-chr path/to/weights/chr \
  --out results/trait
```

The input summary statistics should contain:

```text
SNP A1 A2 N Z
```

If a file contains `BETA` and `SE`, compute `Z = BETA / SE` and write the LDSC
columns before running `ldsc.py`.

## SV sLDSC Reference

The local SV reference is expected at:

```text
ldsc_mixer_sv_reference/
```

This directory was generated with:

```text
make_ldsc_annot.py
make_ldsc.sh
```

`make_ldsc_annot.py` creates LDSC-compatible annotation, partitioned LD-score,
and regression-weight files from the same MiXeR-SV LD matrix and annotation
matrix used by MiXeR-SV. With an SV SNP list, it writes two region annotations:

```text
outside_region
within_region
```

`make_ldsc.sh` is the wrapper used to run `make_ldsc_annot.py` and write the
`ldsc_mixer_sv_reference/` directory. Edit its `LD_DIR`, `ANNOT`, and
`SNP_FILE` variables when rebuilding the reference on a different machine or
from a different MiXeR-SV input directory.

It contains:

```text
ldsc_mixer_sv_reference/ldscore/chr*.l2.ldscore.gz
ldsc_mixer_sv_reference/ldscore/chr*.l2.M
ldsc_mixer_sv_reference/ldscore/chr*.l2.M_5_50
ldsc_mixer_sv_reference/weights/chr*.l2.ldscore.gz
```

The partitioned LD score files include two annotation columns:

```text
outside_regionL2
within_regionL2
```

Thus `--h2` with this reference runs partitioned LDSC/sLDSC and estimates
heritability enrichment for non-SV regions versus SV regions.

Example:

```bash
python ldsc.py \
  --h2 /cluster/projects/nn9114k/datngu/database/alkes_group_data/sumstats_107/PASS.Height.Yengo2022.sumstats.gz \
  --ref-ld-chr ldsc_mixer_sv_reference/ldscore/chr \
  --w-ld-chr ldsc_mixer_sv_reference/weights/chr \
  --out results/height_yengo2022_sv_mixer_test/PASS.Height.Yengo2022.sv_mixer_eur
```

Example Height result:

```text
Total Observed scale h2: 0.6277 (0.0356)
Categories: outside_regionL2_0 within_regionL2_0
Proportion of SNPs: 0.9947 0.0053
Proportion of h2g: 0.9933 0.0067
Enrichment: 0.9986 1.27
```

## Saga Slurm Workflows

Run the Height smoke test:

```bash
sbatch run_height_yengo2022_sv_mixer_test.sh
```

Run sLDSC for all summary statistics in:

```text
/cluster/projects/nn9114k/datngu/database/alkes_group_data/sumstats_107/*sumstats.gz
```

with:

```bash
sbatch run_sldsc_sumstats_107_array.sh
```

The array script uses:

```text
#SBATCH --account=nn9114k
#SBATCH --time=6:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-106%30
```

The manifest is:

```text
manifests/sumstats_107.txt
```

Outputs are written to:

```text
results/sldsc_sumstats_107/
logs/
```

`FIXED.Schizophrenia.Trubetskoy2022.sumstats.gz` uses `ID A1 A2 N BETA SE`
instead of `SNP A1 A2 N Z`. For that file, the completed workflow converts
`ID -> SNP` and `Z = BETA / SE`, then runs:

```bash
sbatch run_sldsc_schizophrenia_converted.sh
```

## Compare sLDSC With MiXeR-SV

MiXeR-SV result files are expected at:

```text
/cluster/projects/nn9114k/datngu/projects/mixer_sv/saga_results_107/*_list.txt
```

Collect and compare the two methods:

```bash
MPLCONFIGDIR=/tmp/ofrei_ldsc_mpl python scripts/collect_sldsc_mixer_sv_results.py \
  --manifest manifests/sumstats_107.txt \
  --ldsc-dir results/sldsc_sumstats_107 \
  --mixer-dir /cluster/projects/nn9114k/datngu/projects/mixer_sv/saga_results_107 \
  --out-dir results/sldsc_sumstats_107/comparison
```

Generated files:

```text
results/sldsc_sumstats_107/comparison/sldsc_mixer_sv_comparison.tsv
results/sldsc_sumstats_107/comparison/sldsc_mixer_sv_correlations.tsv
results/sldsc_sumstats_107/comparison/scatter_within_region_enrichment.png
results/sldsc_sumstats_107/comparison/scatter_total_h2.png
```

Current comparison across 106 traits:

```text
metric                     n    pearson_r    spearman_rho
within_region_enrichment  106  0.77796      0.72685
outside_region_enrichment 106  0.77797      0.71222
total_h2                  106  0.92930      0.98205
```

## Repository Changes From Upstream LDSC

This fork includes:

- Python 3 syntax and runtime compatibility.
- pandas API updates for pandas 1.5.
- package-relative imports in `ldscore`.
- modern NumPy least-squares calls with explicit `rcond=None`.
- Slurm scripts for SV sLDSC workflows.
- scripts to compare sLDSC SV enrichment with MiXeR-SV estimates.

## Citation

If you use LDSC, cite the original LDSC papers:

- Bulik-Sullivan et al. LD Score Regression Distinguishes Confounding from
  Polygenicity in Genome-Wide Association Studies. Nature Genetics, 2015.
- Bulik-Sullivan et al. An Atlas of Genetic Correlations across Human Diseases
  and Traits. Nature Genetics, 2015.
- Finucane et al. Partitioning heritability by functional annotation using
  genome-wide association summary statistics. Nature Genetics, 2015.

## License

LDSC is distributed under the GNU GPL v3 license.
