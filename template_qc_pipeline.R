##-----------------------------------------------------------------------
# QC analysis pipeline
##-----------------------------------------------------------------------

# Load packages
suppressPackageStartupMessages(library(Seurat))
suppressPackageStartupMessages(library(SeuratPipe))
suppressPackageStartupMessages(library(dplyr))
set.seed(12345)
here::here()

##-----------------------------------------------------------------------
# Settings
##-----------------------------------------------------------------------

##-------------------
# Input/Output
io <- list()
io$data_dir <- "your_raw_data_dir_containing_10x_samples/"
io$meta_file <- "sample_metadata.csv"
io$out_dir <- "qc/"
io$tenx_dir <- "/outs/"
io$tenx_counts_dir <- "/filtered_feature_bc_matrix/"


##-------------------
# Opts
opts <- list()
# Threshold for nFeatures
opts$nfeat_thresh <- 1000
# Threshold for mitochondrial content
opts$mito_thresh <- 5
# Columns in metadata file to store in Seurat object
opts$meta_colnames <- c("donor", "condition", "pass_qc")
# QCs to plot (present in Seurat meta.data)
opts$qc_to_plot <- c("nFeature_RNA", "nCount_RNA", "percent.mito")
# Should we use Scrublet for doublet detection?
opts$use_scrublet <- FALSE
# Should we use SoupX for ambient RNA removal?
opts$use_soupx <- FALSE
# Filter genes that are expressed in less than `min.cells`.
opts$min.cells <- 10
# Filter cells that express less than `min.features`.
opts$min.features <- 100
# Filename to store QCed object
opts$obj_filename <- "seu_qc"
# If intermediate seu_preqc.rds file exists, should we force re-analysis
# of Scrublet and SoupX, or read intermediate object? For computing time
# efficiency purposes.
opts$force_reanalysis = FALSE


##-------------------
# Read and filter sample metadata file
sample_meta <- read.csv(file = paste0(io$data_dir, io$meta_file)) |>
  dplyr::filter(condition %in% c("ConditionA", "conditionB")) |>
  dplyr::filter(pass_qc == "TRUE") |>
  dplyr::arrange(sample)


##-----------------------------------------------------------------------
# Analysis pipeline
##-----------------------------------------------------------------------
seu <- run_qc_pipeline(
  data_dir = io$data_dir,
  sample_meta = sample_meta,
  sample_meta_filename = NULL,
  nfeat_thresh = opts$nfeat_thresh,
  mito_thresh = opts$mito_thresh,
  meta_colnames = opts$meta_colnames,
  out_dir = io$out_dir,
  qc_to_plot = opts$qc_to_plot,
  use_scrublet = opts$use_scrublet,
  use_soupx = opts$use_soupx,
  tenx_dir = io$tenx_dir,
  tenx_counts_dir = io$tenx_counts_dir,
  obj_filename = opts$obj_filename,
  force_reanalysis = opts$force_reanalysis,
  min.cells = opts$min.cells,
  min.features = opts$min.features)
