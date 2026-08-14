#### load libraries & utility function
library("Seurat")
library("ggplot2")

# source utility functions
# source("workflow/scripts/utils.R")
# snakemake@source("./utils.R") # does not work when loaded as module (https://github.com/snakemake/snakemake/issues/2205)
source("scripts/utils.R")

# Ported from upstream workflow/scripts/mixscape.R (epigen/mixscape_seurat
# v2.0.3). The upstream script reads snakemake@input/output/config; this port
# receives the same values as positional command-line arguments in a fixed
# order:
#   1  input Seurat object (.rds)
#   2  ALL_object.rds output path (<result_path>/mixscape_seurat/<sample>/ALL_object.rds)
#   3  assay
#   4  variable_features_only (0/1)
#   5  grna_split_symbol
#   6  CalcPerturbSig.grna_col
#   7  CalcPerturbSig.gene_col
#   8  CalcPerturbSig.nt_term
#   9  CalcPerturbSig.n_neighbors
#   10 CalcPerturbSig.ndims
#   11 CalcPerturbSig.split_by_col
#   12 RunMixscape.mixscape_fine_mode
#   13 RunMixscape.lfc_th
#   14 RunMixscape.split_by_col
#   15 RunMixscape.min_de_genes
#   16 RunMixscape.min_cells
#   17 RunMixscape.prtb_type
#   18 MixscapeLDA.npcs
args <- commandArgs(TRUE)
filtered_object_path <- args[1]
mixscape_object_path <- args[2]
sample_dir <- dirname(mixscape_object_path)

# parameters
assay <- args[3]
variable_features_only <- as.integer(args[4])
grna_split_symbol <- args[5]
calcPerturbSig_params <- list(
    grna_col = args[6],
    gene_col = args[7],
    nt_term = args[8],
    n_neighbors = as.integer(args[9]),
    ndims = as.integer(args[10]),
    split_by_col = args[11]
)
runMixscape_params <- list(
    mixscape_fine_mode = args[12],
    lfc_th = as.numeric(args[13]),
    split_by_col = args[14],
    min_de_genes = as.integer(args[15]),
    min_cells = as.integer(args[16]),
    prtb_type = args[17]
)
npcs <- max(calcPerturbSig_params[["ndims"]], as.integer(args[18]))

# outputs (same names as upstream; derived from the sample result dir)
stat_plots_path <- file.path(sample_dir, "plots", "stats")
prtb_data_path <- file.path(sample_dir, "ALL_PRTB_data.csv")
mixscape_stats_path <- file.path(sample_dir, "mixscape_stats.csv")

if (calcPerturbSig_params[["split_by_col"]]==''){
    perturbSig_split_by <- NULL
}else{
    perturbSig_split_by <- calcPerturbSig_params[["split_by_col"]]
}

if (runMixscape_params[["split_by_col"]]==''){
    mixscape_split_by <- NULL
}else{
    mixscape_split_by <- runMixscape_params[["split_by_col"]]
}

# make directories if not exist
if (!dir.exists(stat_plots_path)){
        dir.create(stat_plots_path, recursive = TRUE)
}

### load filtered data
data <- readRDS(file = file.path(filtered_object_path))
DefaultAssay(object = data) <- assay

### Prepare assay for dimensionality reduction: Normalize data, find variable features and scale data.

# check if data is normalized, if not normalize the data
if (all(GetAssayData(object = data, slot = "counts", assay = assay) == GetAssayData(object = data, slot = "data", assay = assay))){
    print("Normalizing data")
    data <- NormalizeData(object = data)
}

# set features for the analysis
if (variable_features_only==1){
    # check if data has variable features, if not determine them
    if (length(VariableFeatures(object = data))==0){
        data <- FindVariableFeatures(object = data)
    }
    features <- VariableFeatures(object = data)
}else{
    features <- rownames(data)
}

# check if data is scaled, if not scale data
if (dim(GetAssayData(object = data, slot = "scale.data", assay = assay))[1]==0){
    print("Scaling data")
    data <- ScaleData(object = data, features = features)
}

# Run Principle Component Analysis (PCA) to reduce the dimensionality of the data.
data <- RunPCA(object = data, features = features, npcs = npcs)

### Mixscape Analysis

# Calculate perturbation signature (PRTB).
# https://satijalab.org/seurat/reference/calcperturbsig
print("Calculating perturbation signature (PRTB)")
data <- CalcPerturbSig(
  object = data,
  assay = assay,
  features = features,
  slot = "data",
  gd.class = calcPerturbSig_params[["gene_col"]],
  nt.cell.class = calcPerturbSig_params[["nt_term"]],
  split.by = perturbSig_split_by,
  num.neighbors = calcPerturbSig_params[["n_neighbors"]],
  reduction = "pca",
  ndims = calcPerturbSig_params[["ndims"]],
  new.assay.name = "PRTB",
  verbose = TRUE
)

# Run mixscape.
# https://satijalab.org/seurat/reference/runmixscape
print("Running Mixscape")
data <- RunMixscape(
  object = data,
  assay = "PRTB",
  slot = "scale.data",
  labels = calcPerturbSig_params[["gene_col"]],
  nt.class.name = calcPerturbSig_params[["nt_term"]],
  new.class.name = "mixscape_class",
  min.de.genes = runMixscape_params[["min_de_genes"]],
  min.cells = runMixscape_params[["min_cells"]],
  de.assay = assay,
  logfc.threshold = runMixscape_params[["lfc_th"]],
  iter.num = 10,
  verbose = TRUE,
  split.by = mixscape_split_by,
  fine.mode = as.logical(runMixscape_params[["mixscape_fine_mode"]]),
  fine.mode.labels = calcPerturbSig_params[["grna_col"]],
  prtb.type = runMixscape_params[["prtb_type"]]
)


### Plot Mixscape Statistics
print("Plot Mixscape statistics")

# Calculate percentage of KO cells for all target gene classes.
stat_table <- table(data$mixscape_class.global, unlist(data[[calcPerturbSig_params[["grna_col"]]]]))
df <- prop.table(stat_table,2)

df2 <- reshape2::melt(df)
df2$Var2 <- as.character(df2$Var2)
test <- df2[which(df2$Var1 == "KO"),]
test <- test[order(test$value, decreasing = T),]
new.levels <- test$Var2
df2$Var2 <- factor(df2$Var2, levels = new.levels )
df2$Var1 <- factor(df2$Var1, levels = c(calcPerturbSig_params[["nt_term"]], "NP", "KO"))
df2$gene <- sapply(as.character(df2$Var2), function(x) head(strsplit(x, split = grna_split_symbol)[[1]],1))
df2$guide_number <- sapply(as.character(df2$Var2), function(x) tail(strsplit(x, split = grna_split_symbol)[[1]],1))
# remove NT
df3 <- df2[-c(which(df2$gene == calcPerturbSig_params[["nt_term"]])),]
df3$Var1 <- factor(df3$Var1, levels = c("NP", "KO"))
df3 <- df3[!is.na(df3$Var1),]

for (gene in unique(df3$gene)){
   tmp_plot <- ggplot(df3[df3$gene==gene,], aes(x = guide_number, y = value*100, fill= Var1)) +
    geom_bar(stat= "identity") +
    theme_classic()+
    scale_fill_manual(values = c("grey79","coral1")) +
    ylab("% of cells") +
    xlab("sgRNA") +
    ggtitle(gene) +
    theme(axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 10),
          strip.text = element_text(size = 10, face = "bold"),
          legend.title = element_text(size = 10),
          legend.text = element_text(size = 10),
         plot.title = element_text(hjust = 0.5)) +
    labs(fill = "mixscape class")

    ggsave_new(filename = gene,
           results_path = stat_plots_path,
           plot = tmp_plot,
           width = 4,
           height = 4)
}

### save results
print("Saving results")
# save seurat object and metadata
save_seurat_object(seurat_obj=data,
                   result_dir=dirname(mixscape_object_path),
                   prefix="ALL_")

# refromat & save mixscape statistics
stat_table <- reshape(as.data.frame(t(stat_table)), timevar = "Var2", idvar = "Var1", direction = "wide")
colnames(stat_table) <- gsub("Freq.", "", colnames(stat_table))
colnames(stat_table)[1] <- "gRNA"
fwrite(as.data.frame(stat_table), file=file.path(mixscape_stats_path), row.names=FALSE)

# save matrix of PRTB values; slots counts and scale.data are emtpy -> only save data slot
fwrite(as.data.frame(GetAssayData(object = data, slot = "data", assay = "PRTB")), file=file.path(prtb_data_path), row.names=TRUE)
