#!/usr/bin/env Rscript
# Generate the tiny Seurat fixture objects used by the oxo-flow mixscape port.
# The fixtures are real, minimal Seurat objects of the same format as the
# upstream workflow's input (a processed per-sample Seurat object with
# gRNAcall/KOcall metadata columns), sized so the dry-run is deterministic
# on CI without downloading real perturbation screens.
#
# Usage: Rscript make_fixtures.R <outdir>   (default: test/fixtures/data)
suppressMessages(library(Seurat))

outdir <- commandArgs(TRUE)[1]
if (is.na(outdir)) outdir <- "test/fixtures/data"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
set.seed(42)

make_fixture <- function(sample_name, outdir) {
    n_genes <- 60
    n_cells <- 30
    counts <- matrix(rpois(n_genes * n_cells, lambda = 1),
                     nrow = n_genes, ncol = n_cells)
    rownames(counts) <- sprintf("GENE%02d", seq_len(n_genes))
    colnames(counts) <- sprintf("%s_CELL%02d", sample_name, seq_len(n_cells))
    obj <- CreateSeuratObject(counts = counts, project = sample_name)
    # metadata columns the Mixscape workflow reads (gRNAcall / KOcall)
    guides <- c(sprintf("STAT1-%d", 1:3), "NonTargeting")
    obj$gRNAcall <- sample(guides, n_cells, replace = TRUE)
    obj$KOcall <- ifelse(obj$gRNAcall == "NonTargeting", "NonTargeting", "STAT1")
    saveRDS(obj, file.path(outdir, paste0(sample_name, ".rds")))
}

for (s in c("S1", "S2")) make_fixture(s, outdir)
cat("wrote:", list.files(outdir, pattern = "\\.rds$"), "\n")
