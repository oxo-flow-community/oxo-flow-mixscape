# oxo-flow mixscape — Pooled CRISPR Perturbation Analysis (Seurat Mixscape)

[![CI](https://github.com/oxo-flow-community/oxo-flow-mixscape/actions/workflows/ci.yml/badge.svg)](https://github.com/oxo-flow-community/oxo-flow-mixscape/actions/workflows/ci.yml)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](LICENSE)

A [Snakemake 8](https://snakemake.readthedocs.io/en/stable/) workflow for
performing perturbation analyses of pooled (multimodal) CRISPR screens with
scRNA-seq read-out (scCRISPR-seq, CROP-seq, Perturb-seq), powered by the R
package [Seurat's](https://satijalab.org/seurat/index.html) method
[Mixscape](https://satijalab.org/seurat/articles/mixscape_vignette.html).

Per input sample (a processed Seurat object), this port calculates local
perturbation signatures (`CalcPerturbSig`), classifies perturbed vs.
non-perturbed cells (`RunMixscape`), analyzes perturbation responses with
Linear Discriminant Analysis (`MixscapeLDA`) and UMAP, produces the
visualization suite (classification statistics, perturbation-score density
plots, posterior-probability violin plots, optional antibody-expression
violin plots), and exports the reproducibility artifacts (conda envs, runtime
config, annotation file).

## Source

Ported from **[epigen/mixscape_seurat](https://github.com/epigen/mixscape_seurat)**,
version `v2.0.3` (MIT). This port is maintained independently and **may lag
the upstream** — check the `v2.0.3` above and the fidelity table below for
the exact ported state. Created 2026-08-15.

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

## Quickstart

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

## Requirements

- **oxo-flow ≥ 0.11.0** — install the prebuilt binary:

```bash
curl -fL -o oxo-flow.tar.gz \
  https://github.com/Traitome/oxo-flow/releases/download/v0.11.0/oxo-flow-v0.11.0-x86_64-unknown-linux-gnu.tar.gz
tar xzf oxo-flow.tar.gz
sudo mv oxo-flow /usr/local/bin/
```

- Conda users may alternatively `conda install -c bioconda oxo-flow-cli`
  (note: the bioconda package currently lags the release binary at 0.10.2 —
  some 0.11.0 format features may not validate).
- Conda at runtime for the environments declared in `main.oxoflow`
  (`envs/seurat_mixscape.yaml`, `envs/seurat_lda.yaml`,
  `envs/config_export.yaml`).

## License

Apache-2.0. Copyright (c) 2026 oxo-flow-community. Upstream attribution in
[NOTICE.md](NOTICE.md).

## Community

https://oxo-flow-community.github.io/
