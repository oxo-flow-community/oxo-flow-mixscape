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
    n_cells <- 60  # >= max(npcs=40, ...) — RunPCA/irlba needs nv < min(nrow, ncol) (live: 'max(nu, nv) must be strictly less than min(nrow(A), ncol(A))' at 30 cells)
    # Gene-varying Poisson means (0.5..20) + two perturbation-marker genes
    # with clearly elevated counts: near-uniform rpois(lambda=1) counts made
    # SCTransform keep zero-variance features only, and RunPCA's irlba died
    # with 'max(nu, nv) must be positive' (live).
    gene_means <- runif(n_genes, 0.5, 20)
    counts <- matrix(rpois(n_genes * n_cells, lambda = rep(gene_means, n_cells)),
                     nrow = n_genes, ncol = n_cells)
    counts[1, ] <- counts[1, ] + rpois(n_cells, lambda = 15)   # marker 1
    counts[2, ] <- counts[2, ] + rpois(n_cells, lambda = 15)   # marker 2
    rownames(counts) <- sprintf("GENE%02d", seq_len(n_genes))
    colnames(counts) <- sprintf("%s_CELL%02d", sample_name, seq_len(n_cells))
    # metadata columns the Mixscape workflow reads (gRNAcall / KOcall)
    guides <- c(sprintf("STAT1-%d", 1:3), "NonTargeting")
    # >=35 NonTargeting cells: CalcPerturbSig's nn2 needs more NT cells
    # than n_neighbors=30 (live: 'Cannot find more nearest neighbours
    # than there are points' when the random draw left ~15 NT cells).
    gRNAcall <- c(rep("NonTargeting", 35),
                   sample(guides[1:3], n_cells - 35, replace = TRUE))
    # Perturbation signature: STAT1 cells overexpress a 5-gene module so
    # CalcPerturbSig has real signal and RunMixscape classifies some
    # cells as KO (live: the plain fixture produced zero KO cells and
    # lda.R's subset(idents = "KO") failed).
    stat1_cells <- which(gRNAcall != "NonTargeting")
    counts[3:12, stat1_cells] <- counts[3:12, stat1_cells] +
        matrix(rpois(10 * length(stat1_cells), lambda = 200),
               nrow = 10, ncol = length(stat1_cells))
    obj <- CreateSeuratObject(counts = counts, project = sample_name)
    obj$gRNAcall <- gRNAcall
    obj$KOcall <- ifelse(gRNAcall == "NonTargeting", "NonTargeting", "STAT1")
    # Upstream (nf-core scrnaseq) feeds Mixscape SCTransform-normalized
    # objects (DefaultAssay = "SCT"); the fixture mirrors that contract.
    obj <- SCTransform(obj, verbose = FALSE)
    saveRDS(obj, file.path(outdir, paste0(sample_name, ".rds")))
}

for (s in c("S1", "S2")) make_fixture(s, outdir)
cat("wrote:", list.files(outdir, pattern = "\\.rds$"), "\n")
