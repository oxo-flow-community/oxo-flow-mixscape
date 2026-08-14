#!/usr/bin/env python3
"""Port of the upstream `config_export` rule (epigen/mixscape_seurat
v2.0.3). The upstream rule is a Snakemake `run:` block that serializes the
runtime config dict as YAML (yaml.dump(config)). oxo-flow renders each
runtime config value into this script's command line, so the exported file
reflects the exact configuration used for the run.

Usage: export_config.py <22 config values in upstream order> <output yaml>
"""
import sys
import yaml

args = sys.argv[1:]
out_path = args[-1]
(
    mem, threads, annotation, result_path, project_name,
    grna_split_symbol, assay, variable_features_only,
    grna_col, gene_col, nt_term, n_neighbors, ndims, cps_split_by_col,
    fine_mode, lfc_th, mixscape_split_by_col, min_de_genes, min_cells,
    prtb_type, lda_npcs, antibody_capture,
) = args[:-1]

# Same structure as the upstream config/config.yaml
config = {
    "mem": mem,
    "threads": threads,
    "annotation": annotation,
    "result_path": result_path,
    "project_name": project_name,
    "grna_split_symbol": grna_split_symbol,
    "assay": assay,
    "variable_features_only": variable_features_only,
    "CalcPerturbSig": {
        "grna_col": grna_col,
        "gene_col": gene_col,
        "nt_term": nt_term,
        "n_neighbors": n_neighbors,
        "ndims": ndims,
        "split_by_col": cps_split_by_col,
    },
    "RunMixscape": {
        "mixscape_fine_mode": fine_mode,
        "lfc_th": lfc_th,
        "split_by_col": mixscape_split_by_col,
        "min_de_genes": min_de_genes,
        "min_cells": min_cells,
        "prtb_type": prtb_type,
    },
    "MixscapeLDA": {
        "npcs": lda_npcs,
    },
    "Antibody_Capture": antibody_capture,
}

with open(out_path, "w") as outfile:
    yaml.dump(config, outfile)
