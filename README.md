# oxo-flow-mixscape — Pooled CRISPR perturbation analysis with Seurat Mixscape

[![CI](https://github.com/oxo-flow-community/oxo-flow-mixscape/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-mixscape/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

> ★ Verified · ⇄ Official port of [`epigen/mixscape_seurat`](https://github.com/epigen/mixscape_seurat) @ `v2.0.3` — same tools, same versions, same commands. Part of the [oxo-flow-community catalog](https://oxo-flow-community.github.io/).

This workflow runs perturbation analyses of pooled (multimodal) CRISPR screens
with scRNA-seq read-out (scCRISPR-seq, CROP-seq, Perturb-seq), powered by the
R package [Seurat's](https://satijalab.org/seurat/index.html) method
[Mixscape](https://satijalab.org/seurat/articles/mixscape_vignette.html).

Starting from one processed Seurat object per sample, it calculates local
perturbation signatures (`CalcPerturbSig`), classifies perturbed vs.
non-perturbed cells (`RunMixscape`), analyzes perturbation responses with
Linear Discriminant Analysis (`MixscapeLDA`) and UMAP, produces the
visualization suite (classification statistics, perturbation-score density
plots, posterior-probability violin plots, optional antibody-expression
violin plots), and exports the reproducibility artifacts (conda envs, runtime
config, annotation file).

## Installation

### 1. Install oxo-flow

Requires **oxo-flow >= 0.12.0**. The recommended way is the prebuilt release
binary (Linux x86_64):

```bash
curl -fL -o oxo-flow.tar.gz https://github.com/Traitome/oxo-flow/releases/latest/download/oxo-flow-latest-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz && sudo mv oxo-flow /usr/local/bin/
```

Alternatively, conda users can install the CLI with
`conda install -c bioconda oxo-flow-cli` — note that the conda package may lag
behind the releases; other platform binaries are available on the
[releases page](https://github.com/Traitome/oxo-flow/releases).

### 2. Get this workflow

```bash
git clone https://github.com/oxo-flow-community/oxo-flow-mixscape.git
```

### 3. Requirements

**Input data (you provide):**

- One processed Seurat object per sample, as `<data_dir>/<sample>.rds` (e.g.
  `data_dir=/path/to/your/data` containing `S1.rds`, `S2.rds`, …). The
  workflow expects already normalized/integrated objects — QC, normalization
  and integration run upstream of this workflow (see the Fidelity section).
- An annotation CSV (`annotation=...`) with at least `name` and `data`
  columns mapping each sample name to its object path (default
  `test/fixtures/annotation.csv`).
- Optional: a 10X `Antibody_Capture` assay named `AB` if you want the
  antibody-expression violin plots (`antibody_capture` config; set to `""`
  to disable).

**Compute:**

- The `mixscape`, `lda` and `visualize` rules each need up to 8 CPUs and
  32000 MB (32 GB) of memory; the export rules need 1 CPU and 1 GB. Use
  `oxo-flow run -j N` to control parallelism.

**Tool delivery:**

- Tools run in conda environments with pinned versions declared in
  `main.oxoflow` — `envs/seurat_mixscape.yaml` and `envs/seurat_lda.yaml`
  (r-seurat 4.4.0, r-seuratobject 4.1.4, r-matrix, r-irlba, r-mixtools,
  r-ggplot2, r-scales, r-patchwork, r-data.table) and `envs/config_export.yaml`
  (pyyaml 6.0.1). Conda (or mamba) must be available at runtime; the
  `env_export_*` rules additionally require a `conda` binary on PATH.

## Usage

```bash
# 1. install oxo-flow (see Requirements)
# 2. prepare data: one processed Seurat .rds per sample, named
#    <data_dir>/<sample>.rds  (see test/fixtures/data/ and annotation.csv)
# 3. preview the plan (defaults point at the fixtures)
oxo-flow dry-run main.oxoflow
# 4. run on your own data
oxo-flow run main.oxoflow -j 8 data_dir=/path/to/your/data annotation=/path/to/your_annotation.csv
# 5. run a subset
oxo-flow run main.oxoflow -t visualize --samples first:2
```

All upstream knobs are exposed as `[config]` keys in `main.oxoflow`
(resources, annotation, data_dir, result_path, project_name, the
`CalcPerturbSig` / `RunMixscape` / `MixscapeLDA` parameters and
`antibody_capture`) and can be overridden on the command line as `key=value`.
Results land under `{config.result_path}/mixscape_seurat/` with the upstream
file names; the `config_export` rule writes the effective runtime
configuration (including CLI overrides) to
`configs/<project_name>_config.yaml`.

## Source

Ported from **[epigen/mixscape_seurat](https://github.com/epigen/mixscape_seurat)** @
`v2.0.3` (MIT, sha `bcf72d5769e2fc7246ba9813892fb438288b6160`). Created
2026-08-15; this workflow may lag behind upstream releases. Full upstream
attribution in [NOTICE.md](NOTICE.md).

## Fidelity

Default-parameters main execution path only. The upstream annotation CSV maps
each sample name to the path of its processed Seurat object; the port reads
the same per-sample `.rds` inputs from `{config.data_dir}/{sample}.rds` and
writes all results under `{config.result_path}/mixscape_seurat/` with the
upstream file names.

| Upstream process/rule | oxo-flow rule | Tool (version) | Notes |
|---|---|---|---|
| `mixscape` | `mixscape` | Seurat 4.4.0 (r-seurat) | identical R script logic (CalcPerturbSig, RunMixscape, stats plots, ALL_* outputs); `snakemake@` I/O/config access replaced with CLI args. Threads 8 = upstream `8 * threads` (threads=1), mem 32000 MB. |
| `lda` | `lda` | Seurat 4.4.0 (r-seurat) | identical R script logic (MixscapeLDA, RunUMAP, FILTERED_*/LDA_data outputs). |
| `visualize` | `visualize` | Seurat 4.4.0 (r-seurat) | identical R script logic (PerturbScore, PosteriorProbability, optional `{Antibody_Capture}_expression` violin plots). The antibody-expression plot dir is produced by the script (same guard as upstream) but is not a declared rule output — oxo-flow cannot declare conditionally-enabled outputs. |
| `env_export` | `env_export_mixscape`, `env_export_lda` | conda (user's install) | upstream's single rule with an `{env}` wildcard split into two explicit rules (environments cannot be wildcarded); `conda env export > {output}` verbatim. `conda` itself is unpinned, exactly as upstream. |
| `config_export` | `config_export` | pyyaml 6.0.1 | upstream `run:` block (`yaml.dump(config)`) ported to `scripts/export_config.py`; runtime config values passed as CLI args. pyyaml pinned at port time (2026-08-15) — upstream relied on the unpinned Snakemake runtime env. |
| `annot_export` | `annot_export` | — | `cp {input} {output}` verbatim. |
| `all` (target) | — | — | not ported: oxo-flow's target is implicit (all rules are targets). |
| demultiplexing | — | — | not ported: upstream module (scrnaseq_processing_seurat), outside this repo. |
| scdna | — | — | not ported: upstream module, outside this repo. |
| normalization (QC/normalization/integration) | — | — | not ported: performed upstream of this workflow by the MrBiomics recipe modules; the input is a processed Seurat object. |
| differential_test (perturbation DE) | — | — | not ported: the per-gene DE runs inside `Seurat::RunMixscape` (min.de.genes / logfc.threshold), not as a separate rule. |

Other deviations: (1) sample input paths come from the `{config.data_dir}`
convention instead of per-row CSV paths — the annotation CSV is retained as
the reproducibility artifact (copied by `annot_export`); (2) the upstream
nested config keys (`CalcPerturbSig.*`, `RunMixscape.*`, `MixscapeLDA.npcs`,
`Antibody_Capture`) are flattened in `[config]` — values and defaults are
identical; (3) `test/fixtures/*.rds` are tiny genuine Seurat objects generated
with Seurat 5.4.0 (local toolchain) for dry-run validation only — upstream
pins r-seurat 4.4.0; (4) the `snakemake@` object access in the R scripts is
replaced with positional CLI args (the ported scripts document the arg
order); (5) a commented-out draft plotting block in upstream `mixscape.R`
was dropped.

## Test

```bash
bash test/run.sh   # validate + lint + dry-run against the fixtures
```

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. The upstream
`epigen/mixscape_seurat` workflow is MIT — see [NOTICE.md](NOTICE.md) and
[LICENSE.upstream](LICENSE.upstream).
