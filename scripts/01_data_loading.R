# =======================================================================
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
# Date: 2025-11-03
#
# Description:
# Loads raw 10X data and metadata, creates individual Seurat objects,
# assigns conditions, and performs initial QC (nFeature/nCount/percent.mt).
#
# Output:
# - A list of Seurat objects saved as RDS (results/qc_filtered/)
# =======================================================================

# ---------------------------
# Load Required Packages
# ---------------------------
suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(Matrix)
  library(ggplot2)
  library(patchwork)
  library(scDblFinder)
  library(SingleCellExperiment)
})

set.seed(42)
options(stringsAsFactors = FALSE)

# ---------------------------
# Define Input Paths
# ---------------------------
# IMPORTANT:
# These paths will be provided by the user when running locally.
# Do NOT hardcode system-specific paths for GitHub.

data_dir     <- "data/raw/"             # expect raw 10X folders here
metadata_csv <- "data/metadata.csv"     # sample metadata provided by user
output_dir   <- "results/qc_filtered/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------
# Step 1 — Load Metadata
# ---------------------------
metadata <- read.csv(metadata_csv, stringsAsFactors = FALSE)

metadata$Sample_ID    <- trimws(metadata$Sample_ID)
metadata$Sample_Title <- trimws(metadata$Sample_Title)
metadata$Condition    <- trimws(metadata$Condition)

print(head(metadata))
print(table(metadata$Condition))   # Should show: NND=6, nAD=6, iAD=13

# ---------------------------
# Step 2 — Verify Folders Match Metadata
# ---------------------------
sample_folders <- list.dirs(path = data_dir, full.names = TRUE, recursive = FALSE)
folder_names   <- basename(sample_folders)

matched <- intersect(metadata$Sample_ID, folder_names)
missing <- setdiff(metadata$Sample_ID, folder_names)

cat("✓ Matched samples:", length(matched), "\n")
if (length(missing) > 0) {
  cat("✗ Missing folders:", paste(missing, collapse = ", "), "\n")
}

# ---------------------------
# Step 3 — Load Raw 10X Samples
# ---------------------------
seurat_list <- list()

for (i in 1:nrow(metadata)) {
  
  sample_id  <- metadata$Sample_ID[i]
  sample_dir <- file.path(data_dir, sample_id)
  
  if (!dir.exists(sample_dir)) {
    warning("Folder not found for: ", sample_id)
    next
  }
  
  message("Loading sample: ", sample_id)
  
  counts <- Read10X(sample_dir)
  
  seu <- CreateSeuratObject(
    counts = counts,
    project = sample_id,
    min.features = 200
  )
  
  seu$Condition     <- metadata$Condition[i]
  seu$Sample_Title  <- metadata$Sample_Title[i]
  
  seurat_list[[sample_id]] <- seu
}

cat("Total samples loaded:", length(seurat_list), "\n")

# ---------------------------
# Step 4 — Calculate QC Metrics
# ---------------------------
for (sample_name in names(seurat_list)) {
  seu <- seurat_list[[sample_name]]
  
  seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")
  
  # Save QC summary per sample (optional)
  message("QC Summary for: ", sample_name)
  print(summary(seu$nFeature_RNA))
  print(summary(seu$nCount_RNA))
  print(summary(seu$percent.mt))
  
  seurat_list[[sample_name]] <- seu
}

# ---------------------------
# Step 5 — Save Output
# ---------------------------
saveRDS(seurat_list, file = file.path(output_dir, "seurat_list_raw_qc.rds"))

cat("\n✓ QC processing completed. Saved to:", output_dir, "\n")
