# =======================================================================
# 06_preintegration_clustering.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Performs clustering BEFORE Harmony integration. This step is used to
#   visualize batch effects and to generate initial PCA/UMAP/HVG diagnostics.
#
# Input:
#   - results/normalized/seurat_list_normalized.rds
#
# Output:
#   - preintegration UMAPs (clusters, samples, conditions)
#   - HVG scatter plot
#   - Elbow plot for 50 PCs
#   Saved under: results/clustering_pre_integration/
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

set.seed(42)

# -----------------------------------------
# Define IO paths
# -----------------------------------------
input_rds  <- "results/normalized/seurat_list_normalized.rds"
output_dir <- "results/clustering_pre_integration/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------
# Load normalized samples
# -----------------------------------------
seurat_normalized <- readRDS(input_rds)

cat("Loaded", length(seurat_normalized), "normalized samples.\n")

# ============================================================================
# STEP 1 — Select top 5000 variable features
# ============================================================================
cat("Selecting top 5000 integration features...\n")

features <- SelectIntegrationFeatures(
  object.list = seurat_normalized,
  nfeatures = 5000
)

cat("Selected", length(features), "variable features.\n")

# ============================================================================
# STEP 2 — Merge samples (NO integration)
# ============================================================================
cat("Merging samples (pre-integration)...\n")

combined_raw <- merge(
  seurat_normalized[[1]],
  y = seurat_normalized[-1],
  merge.data = TRUE
)

cat("Merged cell count:", ncol(combined_raw), "\n")

# ============================================================================
# STEP 3 — PCA (50 PCs)
# ============================================================================
cat("Running PCA...\n")

combined_raw <- RunPCA(
  combined_raw,
  features = features,
  npcs = 50,
  verbose = FALSE
)

# ============================================================================
# STEP 4 — UMAP (Dims 1–30)
# ============================================================================
cat("Running UMAP...\n")

combined_raw <- RunUMAP(
  combined_raw,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

# ============================================================================
# STEP 5 — Clustering
# ============================================================================
cat("Finding neighbors (dims 1–30)...\n")

combined_raw <- FindNeighbors(
  combined_raw,
  reduction = "pca",
  dims = 1:30,
  verbose = FALSE
)

combined_raw <- FindClusters(
  combined_raw,
  resolution = 0.5,
  verbose = FALSE
)

cat("\n✓ Clustering complete!\n")

# ============================================================================
# VISUAL 1 — UMAP by cluster
# ============================================================================
p_clusters <- DimPlot(
  combined_raw,
  reduction = "umap",
  label = TRUE,
  label.size = 5,
  pt.size = 0.8
) +
  ggtitle("Clusters BEFORE Integration") +
  theme(plot.title = element_text(size = 16, face = "bold"))

ggsave(
  file.path(output_dir, "clusters_before_integration.jpeg"),
  p_clusters, width = 7, height = 6, dpi = 300
)

# ============================================================================
# VISUAL 2 — UMAP by sample (batch effect)
# ============================================================================
p_sample <- DimPlot(
  combined_raw,
  reduction = "umap",
  group.by = "orig.ident",
  pt.size = 0.6
) +
  ggtitle("Samples BEFORE Integration (Batch Effect Visible)")

ggsave(
  file.path(output_dir, "samples_before_integration.jpeg"),
  p_sample, width = 7, height = 6, dpi = 300
)

# ============================================================================
# VISUAL 3 — UMAP by condition
# ============================================================================
p_condition <- DimPlot(
  combined_raw,
  reduction = "umap",
  group.by = "Condition",
  pt.size = 0.6,
  cols = c("NND" = "#3498db", "nAD" = "#e74c3c", "iAD" = "#2ecc71")
) +
  ggtitle("Conditions BEFORE Integration")

ggsave(
  file.path(output_dir, "conditions_before_integration.jpeg"),
  p_condition, width = 7, height = 6, dpi = 300
)

# ============================================================================
# HVG Scatter Plot (RNA assay)
# ============================================================================
cat("Generating HVG scatter plot...\n")

DefaultAssay(combined_raw) <- "RNA"

combined_raw <- FindVariableFeatures(
  combined_raw,
  selection.method = "vst",
  nfeatures = 2000
)

hv_plot <- VariableFeaturePlot(combined_raw)
top20   <- head(VariableFeatures(combined_raw), 20)

hv_labeled <- LabelPoints(
  plot = hv_plot,
  points = top20,
  repel = TRUE
) + ggtitle("Highly Variable Genes (RNA)")

ggsave(
  file.path(output_dir, "HVG_scatter_plot.jpeg"),
  hv_labeled, width = 7, height = 6, dpi = 300
)

# ============================================================================
# Elbow Plot (50 PCs)
# ============================================================================
cat("Saving elbow plot...\n")

DefaultAssay(combined_raw) <- "SCT"  # restore SCT

temp_pca <- RunPCA(
  combined_raw,
  npcs = 50,
  verbose = FALSE,
  features = VariableFeatures(combined_raw)
)

elbow_plot <- ElbowPlot(temp_pca, ndims = 50) +
  ggtitle("Elbow Plot (50 PCs)") +
  theme_bw()

ggsave(
  file.path(output_dir, "elbow_plot_50PCs.jpeg"),
  elbow_plot, width = 7, height = 6, dpi = 300
)

cat("\n✓ All PRE-INTEGRATION clustering outputs saved successfully.\n")
