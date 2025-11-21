# =======================================================================
# 02_qc_pre_filtering.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
# Description:
#   Performs QC visualization prior to filtering. Generates violin plots,
#   scatter plots, and histograms of key QC metrics.
#
# Input:
#   - seurat_list_raw_qc.rds (output from 01_data_loading.R)
#
# Output:
#   - Violin plots (nFeature, nCount, percent.mt)
#   - QC scatter plots
#   - Histograms
#   Saved under: results/qc_plots/before_filtering/
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

set.seed(42)

# ---------------------------
# Define Paths
# ---------------------------
input_rds  <- "results/qc_filtered/seurat_list_raw_qc.rds"
output_dir <- "results/qc_plots/before_filtering/"

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# ---------------------------
# Load Data
# ---------------------------
seurat_list <- readRDS(input_rds)

# Merge samples for global QC visualization
combined_raw <- merge(seurat_list[[1]], y = seurat_list[-1])

# Add mitochondrial percentage
combined_raw[["percent.mt"]] <- PercentageFeatureSet(combined_raw, pattern = "^MT-")

# --------------------------------
# Violin Plots (nFeature / nCount / percent.mt)
# --------------------------------

clean_theme <- theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(angle = 90, hjust = 1, size = 6),
    legend.position = "none"
  )

vln_nFeature <- VlnPlot(
  combined_raw, features = "nFeature_RNA", group.by = "orig.ident",
  pt.size = 0
) + clean_theme

vln_nCount <- VlnPlot(
  combined_raw, features = "nCount_RNA", group.by = "orig.ident",
  pt.size = 0
) + clean_theme

vln_percent <- VlnPlot(
  combined_raw, features = "percent.mt", group.by = "orig.ident",
  pt.size = 0
) + clean_theme

# Save violin plots
ggsave(file.path(output_dir, "before_qc_vln_nFeature.jpeg"), vln_nFeature, width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "before_qc_vln_nCount.jpeg"),   vln_nCount,   width = 10, height = 4, dpi = 300)
ggsave(file.path(output_dir, "before_qc_vln_percentMT.jpeg"),vln_percent,  width = 10, height = 4, dpi = 300)

# --------------------------------
# QC Scatter Plots
# --------------------------------

qc_df <- combined_raw@meta.data %>%
  mutate(orig.ident = as.factor(orig.ident))

p1 <- ggplot(qc_df, aes(nCount_RNA, nFeature_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw(base_size = 14) +
  labs(title = "nCount_RNA vs nFeature_RNA (Before Filtering)") +
  theme(legend.position = "none")

p2 <- ggplot(qc_df, aes(percent.mt, nFeature_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw(base_size = 14) +
  labs(title = "percent.mt vs nFeature_RNA") +
  theme(legend.position = "none")

p3 <- ggplot(qc_df, aes(percent.mt, nCount_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw(base_size = 14) +
  labs(title = "percent.mt vs nCount_RNA") +
  theme(legend.position = "none")

p4 <- ggplot(qc_df, aes(percent.mt, fill = orig.ident)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  theme_bw(base_size = 14) +
  labs(title = "Histogram of percent.mt") +
  theme(legend.position = "none")

ggsave(file.path(output_dir, "QC_nCount_vs_nFeature.jpeg"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "QC_percentMT_vs_nFeature.jpeg"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "QC_percentMT_vs_nCount.jpeg"), p3, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir, "QC_percentMT_histogram.jpeg"), p4, width = 6, height = 5, dpi = 300)

cat("\n✓ QC BEFORE filtering plots saved to:", output_dir, "\n")
