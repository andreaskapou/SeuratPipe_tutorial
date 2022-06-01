##-----------------------------------------------------------------------
# Independent clustering analysis pipeline
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
set.seed(12345)
here::here()

##-----------------------------------------------------------------------
# Load QCed data
##-----------------------------------------------------------------------
data_dir <- "qc/" # Folder with QC output
obj <- readRDS(paste0(data_dir, "/seu_qc.rds"))
seu <- obj$seu
opts <- obj$opts # Opts keeps information about previous analysis steps for reproducibility
rm(obj)
gc(verbose = FALSE)


##-----------------------------------------------------------------------
# Additional Settings
##-----------------------------------------------------------------------

##-------------------
# Input/Output
io <- list()
io$out_dir <- "cluster/"


##-------------------
# Opts
# Number of PCs, can be a vector: c(30, 50)
opts$npcs <- c(50)
# Top PCs to perform UMAP and clustering, can be a vector: c(30, 40)
opts$ndims <- c(30)
# Clustering resolutions
opts$res <- seq(from = 0.1, to = 0.3, by = 0.1)
# QC information to plot
opts$qc_to_plot <- c("nFeature_RNA", "percent.mito")
# Metadata columns to plot
opts$meta_to_plot <- c("condition")


# Test genes that are detected in a minimum fraction of min.pct cells
opts$min.pct <- 0.25
# Test genes that show, on average, at least X-fold difference
# between the two groups of cells.
opts$logfc.threshold <- 0.5
# Only return positive markers
opts$only.pos <- TRUE
# Maximum number of marker genes to plot for each cluster
opts$topn_genes <- 10
# Should we create feature plots of cluster markers?
opts$plot_cluster_markers = FALSE


# Which PCs to remove prior to running harmony (technical effects)
# Can be a vector e.g. c(1, 2, 4). If NULL, all PCs are used as input.
opts$pcs_to_remove <- NULL
# Maximum cutoff values for plotting continuous features, e.g. gene expression
# Gives better plots where colour scale is not driven by a (few) outlier cells.
# Set to NULL to recoved default Seurat plots
opts$max.cutoff <- "q98"


# Number of highly variable genes to compute
opts$n_hvgs <- 3000
# Set specific seed for reproducibility
opts$seed <- 1
# Discrete colour palette
opts$discrete_col_pal <- SeuratPipe:::discrete_col_pal
# Whether to label the clusters in 'plot_reduction' space.
opts$label <- TRUE
# Sets size of labels.
opts$label.size <- 8
# Adjust point size for plotting.
opts$pt.size <- 1.5
# Figure resolution in ppi
opts$fig.res = 200


##-------------------
# Load modules

# For instance define modules for lineage and cell cycle
lineage <- list(lin_CD4 = c("IL7R", "CCR7"),
                lin_Monocyte = c("CD14", "LYZ", "FCGR3A"),
                lin_NK = c("GNLY", "NKG7"))
# Take a random subset from Seurat's cycle genes
cell_cycle <- list(cc_s = c("PCNA", "RRM1", "POLA1", "USP1"),
                   cc_g2m = c("HMGB2", "MKI67", "TMPO", "CTCF"))

# Group all modules in named list to pass to SeuratPipe functions
opts$modules_group <- list(lineage = lineage,
                           cell_cycle = cell_cycle)


##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------
for (s in names(seu)) {
  obj <- seu[[s]]
  # If sample didn't pass QC
  if (is.null(obj)) { next }
  # If sample contains small number of cells, probabaly ignore it from analysis?
  if (NCOL(obj) < 100) { next }
  # Set analysis directory for each sample
  sample_dir <- paste0(io$out_dir, "/", s, "/")

  seu[[s]] <- run_cluster_pipeline(
    seu_obj = obj,
    out_dir = sample_dir,
    npcs = opts$npcs,
    ndims = opts$ndims,
    res = opts$res,
    modules_group = opts$modules_group,
    metadata_to_plot = opts$meta_to_plot,
    qc_to_plot = opts$qc_to_plot,
    logfc.threshold = opts$logfc.threshold,
    min.pct = opts$min.pct,
    only.pos = opts$only.pos,
    topn_genes = opts$topn_genes,
    pcs_to_remove = opts$pcs_to_remove,
    plot_cluster_markers = opts$plot_cluster_markers,
    max.cutoff = opts$max.cutoff,
    n_hvgs = opts$n_hvgs,
    seed = opts$seed,
    discrete_col_pal = opts$discrete_col_pal,
    label = opts$label,
    label.size = opts$label.size,
    pt.size = opts$pt.size,
    fig.res = opts$fig.res)
}
