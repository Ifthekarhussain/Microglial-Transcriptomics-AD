# =======================================================================
# 07_harmony_integration.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Performs PCA + Harmony batch correction + UMAP + initial clustering.
#   Outputs a fully integrated Seurat object.
#
# Input:
#   - results/normalized/seurat_list_normalized.rds
#
# Output:
#   - results/integrated/combined_integrated.rds
#   - UMAP plots (clusters, samples, conditions)
#
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------
# Define paths
# -----------------------------------------
input_rds   <- "results/normalized/seurat_list_normalized.rds"
output_rds  <- "results/integrated/combined_integrated.rds"
plot_dir    <- "results/umap_post_integration/"

dir.create("results/integrated/", recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------
# Load SCT normalized samples
# -----------------------------------------
seurat_norm <- readRDS(input_rds)

cat("Loaded", length(seurat_norm), "normalized samples.\n")

# -----------------------------------------
# Merge samples (no integration yet)
# -----------------------------------------
combined <- merge(
  seurat_norm[[1]],
  y = seurat_norm[-1],
  merge.data = TRUE
)

DefaultAssay(combined) <- "SCT"
combined <- UpdateSeuratObject(combined)

# Use HVGs from first sample (standard practice)
VariableFeatures(combined) <- VariableFeatures(seurat_norm[[1]])

# -----------------------------------------
# PCA
# -----------------------------------------
cat("Running PCA...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# -----------------------------------------
# HARMONY BATCH CORRECTION
# -----------------------------------------
cat("Running Harmony integration...\n")

combined_integrated <- RunHarmony(
  object = combined,
  dims.use = 1:30,
  group.by.vars = "orig.ident",
  theta = 2,
  plot_convergence = FALSE,
  verbose = FALSE,
  assay.use = "SCT"
)

# -----------------------------------------
# UMAP
# -----------------------------------------
cat("Running UMAP...\n")

combined_integrated <- RunUMAP(
  combined_integrated,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

# -----------------------------------------
# Neighbors + Clustering (default res = 0.5)
# -----------------------------------------
cat("Running graph clustering...\n")

combined_integrated <- FindNeighbors(
  combined_integrated,
  reduction = "harmony",
  dims = 1:30,
  verbose = FALSE
)

combined_integrated <- FindClusters(
  combined_integrated,
  resolution = 0.5,
  verbose = FALSE
)

cat("\n✓ Harmony Integration Complete!\n")
cat("Cells:", ncol(combined_integrated), "\n")
cat("Clusters:", length(unique(combined_integrated$seurat_clusters)), "\n")

# -----------------------------------------
# Save integrated object
# -----------------------------------------
saveRDS(combined_integrated, output_rds)

# -----------------------------------------
# Plots
# -----------------------------------------

# Clusters
p_clusters <- DimPlot(
  combined_integrated,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5,
  label.size = 5
) +
  ggtitle("Clusters AFTER Harmony Integration")

ggsave(file.path(plot_dir, "clusters_after_integration.jpeg"),
       p_clusters, width = 7, height = 6, dpi = 300)

# By sample
p_sample <- DimPlot(
  combined_integrated,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.4
) +
  ggtitle("Samples AFTER Integration (Batch Mixed)")

ggsave(file.path(plot_dir, "samples_after_integration.jpeg"),
       p_sample, width = 7, height = 6, dpi = 300)

# By Condition
p_condition <- DimPlot(
  combined_integrated,
  reduction = "umap",
  group.by = "Condition",
  pt.size = 0.6,
  cols = c("NND" = "#1f78b4", "nAD" = "#33a02c", "iAD" = "#e31a1c")
) +
  ggtitle("Conditions AFTER Integration")

ggsave(file.path(plot_dir, "condition_after_integration.jpeg"),
       p_condition, width = 7, height = 6, dpi = 300)

cat("\n✓ Saved all Harmony integration outputs.\n")
