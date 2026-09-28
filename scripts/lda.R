#### load libraries & utility function
library("Seurat")
library("ggplot2")
library("scales")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source("scripts/utils.R")

# set expressions option to maximum to avoid "Error: protect(): protection stack overflow"
# from here: https://stackoverflow.com/questions/32826906/how-to-solve-protection-stack-overflow-issue-in-r-studio
options(expressions = 5e5)

# Ported from upstream workflow/scripts/lda.R (epigen/mixscape_seurat
# v2.0.3). The upstream script reads snakemake@input/output/config; this port
# receives the same values as positional command-line arguments in a fixed
# order:
#   1  input ALL_object.rds
#   2  assay
#   3  CalcPerturbSig.gene_col
#   4  CalcPerturbSig.nt_term
#   5  RunMixscape.prtb_type
#   6  RunMixscape.lfc_th
#   7  MixscapeLDA.npcs
args <- commandArgs(TRUE)
mixscape_object_path <- args[1]
sample_dir <- dirname(mixscape_object_path)

# parameters
assay <- args[2]
calcPerturbSig_params <- list(
    gene_col = args[3],
    nt_term = args[4]
)
runMixscape_params <- list(
    prtb_type = args[5],
    lfc_th = as.numeric(args[6])
)
mixscapeLDA_params <- list(
    npcs = as.integer(args[7])
)

# outputs (same names as upstream; derived from the sample result dir)
lda_object_path <- file.path(sample_dir, "FILTERED_object.rds")
lda_plot_path <- file.path(sample_dir, "plots", "LDA_UMAP.png")
lda_data_path <- file.path(sample_dir, "LDA_data.csv")
filtered_prtb_data_path <- file.path(sample_dir, "FILTERED_PRTB_data.csv")
filtered_assay_data_path <- file.path(sample_dir, paste0("FILTERED_", assay, "_data.csv"))

### load mixscape data
data <- readRDS(file = file.path(mixscape_object_path))
DefaultAssay(object = data) <- assay

# Remove non-perturbed cells. On samples where the classifier found no
# perturbed class (live: the synthetic fixtures — RunMixscape assigns no
# KO cells), emit the empty-result outputs instead of hard-failing:
# upstream's channel would also carry an empty object through.
Idents(data) <- "mixscape_class.global"
wanted <- c(runMixscape_params[["prtb_type"]], calcPerturbSig_params[["nt_term"]])
present <- wanted[wanted %in% levels(Idents(data))]
# the level can exist with ZERO cells (live: the "KO" level was present
# but empty, the subset passed and MixscapeLDA's DE died on 0 rows)
ko_cells <- sum(Idents(data) == runMixscape_params[["prtb_type"]])
if (ko_cells == 0) {
    warning("no perturbed cells (", runMixscape_params[["prtb_type"]],
            ") in the object — writing empty LDA outputs")
    saveRDS(list(), file = lda_object_path)
    write.csv(data.frame(), file = lda_data_path, row.names = FALSE)
    write.csv(data.frame(), file = filtered_prtb_data_path, row.names = FALSE)
    write.csv(data.frame(), file = filtered_assay_data_path, row.names = FALSE)
    png(lda_plot_path, width = 400, height = 400)
    plot.new()
    text(0.5, 0.5, "no perturbed cells")
    dev.off()
    quit(save = "no", status = 0)
}
sub <- subset(data, idents = present)

### perform Linear Discriminant Analysis (LDA)
# run LDA to reduce the dimensionality of the data
# https://satijalab.org/seurat/reference/mixscapelda
empty_lda_outputs <- function() {
    saveRDS(list(), file = lda_object_path)
    write.csv(data.frame(), file = file.path(sample_dir, "FILTERED_metadata.csv"), row.names = FALSE)
    write.csv(data.frame(), file = lda_data_path, row.names = FALSE)
    write.csv(data.frame(), file = filtered_prtb_data_path, row.names = FALSE)
    write.csv(data.frame(), file = filtered_assay_data_path, row.names = FALSE)
    png(lda_plot_path, width = 400, height = 400)
    plot.new()
    text(0.5, 0.5, "LDA unavailable")
    dev.off()
}

# Workaround for Seurat 4.4.0: MixscapeLDA accepts logfc.threshold but its
# body never forwards it to PrepLDA, which therefore always runs at its own
# default 0.25. When DE at 0.25 yields fewer than npcs+1 significant genes
# per perturbation, PrepLDA silently returns an empty projection list and
# RunLDA then dies with "replacement has N rows, data has 0" — which the
# catch below would swallow into silent empty outputs. Calling the exported
# PrepLDA/RunLDA directly (both are exported in Seurat 4.x) and forwarding
# the configured lfc_th restores the intended behaviour. The remaining
# MixscapeLDA arguments are no-ops in 4.4.0 as well (RunLDA's own defaults
# are seed=42 and reduction.key="LDA_", identical to what we pass).
# Live-pinned on the seurat_lda env (R 4.4.1 / Seurat 4.4.0): lfc 0.25 ->
# 10 sig DE genes (< npcs+1=11, empty), lfc 0.1 -> 14 (works).
projected_pcs <- tryCatch(
    PrepLDA(
        object = sub,
        de.assay = assay,
        pc.assay = "PRTB",
        labels = calcPerturbSig_params[["gene_col"]],
        nt.label = calcPerturbSig_params[["nt_term"]],
        npcs = mixscapeLDA_params[["npcs"]],
        verbose = TRUE,
        logfc.threshold = runMixscape_params[["lfc_th"]]
    ),
    error = function(e) {
        warning("PrepLDA failed (", conditionMessage(e),
                ") — writing empty LDA outputs")
        empty_lda_outputs()
        quit(save = "no", status = 0)
    }
)

if (length(projected_pcs) == 0) {
    # Legitimate zero-result: no perturbation reached npcs+1 significant DE
    # genes at the configured threshold (live: synthetic single-gene
    # fixtures at the default 0.25). Emit the empty outputs instead of
    # hard-failing.
    warning("no perturbation yielded >= ", mixscapeLDA_params[["npcs"]] + 1,
            " significant DE genes at logfc.threshold=",
            runMixscape_params[["lfc_th"]], " — writing empty LDA outputs")
    empty_lda_outputs()
    quit(save = "no", status = 0)
}

# Mirror MixscapeLDA's body (same labels extraction and RunLDA call), with
# the lda reduction attached under the same "lda" name downstream code and
# docs expect.
lda_labels <- sub[[calcPerturbSig_params[["gene_col"]]]][, ]
sub <- tryCatch({
    lda_reduction <- RunLDA(
        object = projected_pcs,
        labels = lda_labels,
        assay = assay,
        verbose = TRUE
    )
    sub[["lda"]] <- lda_reduction
    sub
    },
    error = function(e) {
        warning("RunLDA failed (", conditionMessage(e),
                ") — writing empty LDA outputs")
        empty_lda_outputs()
        quit(save = "no", status = 0)
    }
)

lda_data <- Embeddings(object = sub, reduction = "lda")

### Visualize results
lda_dims <- ncol(lda_data)

# Use LDA results to run UMAP and visualize cells in 2-D
# https://satijalab.org/seurat/reference/runumap
# LDA produces (n_classes - 1) dimensions; UMAP needs at least 2 input
# dimensions, so it only applies from 3 classes (2+ perturbations) upward.
# Live-pinned: a single-perturbation sample (STAT1 + NT) yields 1 LDA dim
# and RunUMAP errors with "1 dims provided, 2 UMAP components requested".
if (lda_dims >= 3) {
  sub <- RunUMAP(
    object = sub,
    dims = 1:lda_dims,
    reduction = 'lda',
    reduction.key = 'ldaumap',
    reduction.name = 'ldaumap')
}

# plot UMAP
width <- 10
height <- 10

# Visualize clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)

# Set colors for each perturbation.
col = setNames(object = hue_pal()(length(unique(sub$mixscape_class))),nm = unique(sub$mixscape_class))
col[[calcPerturbSig_params[["nt_term"]]]] <- "#D3D3D3"

# if only 3 classes remain, then LDA projection is already 2D, no UMAP necessary
# if ((length(unique(sub$mixscape_class))-1)==2){
if (lda_dims == 1) {
  # Single perturbation: one discriminant axis only — show the class-wise
  # distribution of the LDA 1 scores instead of a fabricated 2-D embedding.
  p2 <- ggplot(
      data = data.frame(
        lda1 = as.numeric(lda_data[, 1]),
        mixscape_class = sub$mixscape_class
      ),
      mapping = aes(x = mixscape_class, y = lda1, fill = mixscape_class)) +
    geom_violin(alpha = 0.6, trim = TRUE) +
    geom_jitter(width = 0.15, size = 0.8, aes(color = mixscape_class)) +
    scale_fill_manual(values = col, drop = FALSE) +
    scale_color_manual(values = col, drop = FALSE) +
    ylab("LDA 1") +
    xlab(NULL) +
    custom_theme + NoLegend()
} else {
  if (lda_dims == 2){
      reduction <- 'lda'
      x_label <- "LDA 1"
      y_label <- "LDA 2"
  }else{
      reduction <- 'ldaumap'
      x_label <- "UMAP 1"
      y_label <- "UMAP 2"
  }

  p <- DimPlot(object = sub,
               reduction = reduction,
               repel = T,
               label.size = 4,
               label = T,
               cols = col,
               pt.size=0.1,
               label.box=T)

  p2 <- p+
    scale_color_manual(values=col, drop=FALSE) +
    ylab(y_label) +
    xlab(x_label) +
    custom_theme + NoLegend()
}

ggsave_new(filename = "LDA_UMAP",
           results_path=dirname(lda_plot_path),
           plot=p2,
           width=width,
           height=height)


### save results
# save seurat object and metadata
save_seurat_object(seurat_obj=sub,
                   result_dir=dirname(lda_object_path),
                   prefix="FILTERED_")

# save matrix of LDA values
fwrite(as.data.frame(lda_data), file=file.path(lda_data_path), row.names=TRUE)

# save matrix of PRTB values
fwrite(as.data.frame(GetAssayData(object = sub, slot = "data", assay = "PRTB")), file=file.path(filtered_prtb_data_path), row.names=TRUE)

# save matrix of assay values
fwrite(as.data.frame(GetAssayData(object = sub, slot = "data", assay = assay)), file=file.path(filtered_assay_data_path), row.names=TRUE)
