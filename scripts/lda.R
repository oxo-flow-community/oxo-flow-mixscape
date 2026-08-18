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
if (length(setdiff(wanted, present)) > 0 && length(present) == 1) {
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
sub <- MixscapeLDA(
  object = sub,
  assay = assay,
  ndims.print = 1:5,
  nfeatures.print = 30,
  reduction.key = "LDA_",
  seed = 42,
  pc.assay = "PRTB",
  labels = calcPerturbSig_params[["gene_col"]],
  nt.label = calcPerturbSig_params[["nt_term"]],
  npcs = mixscapeLDA_params[["npcs"]],
  verbose = TRUE,
  logfc.threshold = runMixscape_params[["lfc_th"]]
)

lda_data <- Embeddings(object = sub, reduction = "lda")

### Visualize results

# Use LDA results to run UMAP and visualize cells in 2-D
# https://satijalab.org/seurat/reference/runumap
sub <- RunUMAP(
  object = sub,
  dims = 1:ncol(lda_data),#(length(unique(sub$mixscape_class))-1),
  reduction = 'lda',
  reduction.key = 'ldaumap',
  reduction.name = 'ldaumap')

# plot UMAP
width <- 10
height <- 10


# if only 3 classes remain, then LDA projection is already 2D, no UMAP necessary
# if ((length(unique(sub$mixscape_class))-1)==2){
if (ncol(lda_data)==2){
    reduction <- 'lda'
    x_label <- "LDA 1"
    y_label <- "LDA 2"
}else{
    reduction <- 'ldaumap'
    x_label <- "UMAP 1"
    y_label <- "UMAP 2"
}

# Visualize UMAP clustering results.
Idents(sub) <- "mixscape_class"
sub$mixscape_class <- as.factor(sub$mixscape_class)

# Set colors for each perturbation.
col = setNames(object = hue_pal()(length(unique(sub$mixscape_class))),nm = unique(sub$mixscape_class))
col[[calcPerturbSig_params[["nt_term"]]]] <- "#D3D3D3"

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
