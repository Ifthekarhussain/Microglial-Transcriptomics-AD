# =======================================================================
# 03_qc_filtering.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Applies QC filtering per sample using MAD-based thresholds for
#   nCount_RNA, nFeature_RNA, and percent.mt. Generates BEFORE and AFTER
#   QC violin plots + scatter plots.
#
# Input:
#   - seurat_list_raw_qc.rds (from 01_data_loading.R)
#
# Output:
#   - QC-filtered Seurat list (seurat_list_qc.rds)
#   - BEFORE/AFTER QC plots saved under:
#       results/qc_plots/before_filtering/
#       results/qc_plots/after_filtering/
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# ---------------------------
# Paths
# ---------------------------
input_rds  <- "results/qc_filtered/seurat_list_raw_qc.rds"
output_dir_before <- "results/qc_plots/before_filtering/"
output_dir_after  <- "results/qc_plots/after_filtering/"
output_qc_rds     <- "results/qc_filtered/seurat_list_qc.rds"

dir.create(output_dir_before, recursive = TRUE, showWarnings = FALSE)
dir.create(output_dir_after,  recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load input data
# ---------------------------
seurat_list <- readRDS(input_rds)

cat("Loaded samples:", length(seurat_list), "\n")

# =======================================================================
# STEP 1 — Add percent.mt
# =======================================================================
cat("Adding percent.mt...\n")

for (sid in names(seurat_list)) {
  seurat_list[[sid]][["percent.mt"]] <- PercentageFeatureSet(
    seurat_list[[sid]], pattern = "^MT-"
  )
}

cat("✓ Added percent.mt\n")

# =======================================================================
# STEP 2 — Compute QC thresholds
# =======================================================================

qc_thresholds <- lapply(seurat_list, function(obj) {
  tibble(
    Sample_ID       = obj@project.name,
    median_nCount   = median(obj$nCount_RNA),
    mad_nCount      = mad(obj$nCount_RNA),
    median_nFeature = median(obj$nFeature_RNA),
    mad_nFeature    = mad(obj$nFeature_RNA)
  )
}) %>%
  bind_rows() %>%
  mutate(
    thr_nCount    = pmax(0, median_nCount - 3 * mad_nCount),
    thr_nFeature  = pmax(0, median_nFeature - 2 * mad_nFeature),
    thr_percentMT = 20
  )

cat("✓ QC thresholds computed\n")
print(head(qc_thresholds))

# =======================================================================
# STEP 3 — Apply QC filters
# =======================================================================

seurat_qc <- list()

for (sid in names(seurat_list)) {
  obj <- seurat_list[[sid]]
  thr <- qc_thresholds %>% filter(Sample_ID == sid)

  obj_f <- subset(
    obj,
    subset = nCount_RNA > thr$thr_nCount &
             nFeature_RNA > thr$thr_nFeature &
             percent.mt < thr$thr_percentMT
  )

  seurat_qc[[sid]] <- obj_f
}

cat("✓ QC filtering complete\n")

# =======================================================================
# STEP 4 — Summary Table
# =======================================================================

qc_summary <- tibble(
  Sample_ID      = names(seurat_list),
  Cells_Before   = sapply(seurat_list, ncol),
  Cells_After    = sapply(seurat_qc, ncol),
  Condition      = sapply(seurat_list, \(x) unique(x$Condition)[1])
) %>%
  mutate(
    Cells_Removed    = Cells_Before - Cells_After,
    Percent_Retained = round(100 * Cells_After / Cells_Before, 1)
  )

print(qc_summary)

cat("\nTotal Before QC:", sum(qc_summary$Cells_Before))
cat("\nTotal After QC:",  sum(qc_summary$Cells_After))
cat("\nRetention:", round(100 * sum(qc_summary$Cells_After) /
                         sum(qc_summary$Cells_Before), 1), "%\n\n")

saveRDS(seurat_qc, output_qc_rds)

# =======================================================================
# STEP 5 — BEFORE QC plots
# =======================================================================

combined_before <- merge(seurat_list[[1]], y = seurat_list[-1])
combined_before$orig.ident <- as.factor(combined_before$orig.ident)

# Violin plot: 3 metrics
p_before_vln <- VlnPlot(
  combined_before,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  ncol = 3,
  pt.size = 0.1
) + theme_bw(base_size = 10)

ggsave(
  file.path(output_dir_before, "beforeQC_three_metrics.jpeg"),
  p_before_vln, width = 12, height = 4, dpi = 300
)

# =======================================================================
# STEP 6 — AFTER QC plots
# =======================================================================

combined_after <- merge(seurat_qc[[1]], y = seurat_qc[-1])
combined_after$orig.ident <- as.factor(combined_after$orig.ident)

p_after_vln <- VlnPlot(
  combined_after,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  group.by = "orig.ident",
  ncol = 3,
  pt.size = 0.1
) + theme_bw(base_size = 10)

ggsave(
  file.path(output_dir_after, "afterQC_three_metrics.jpeg"),
  p_after_vln, width = 12, height = 4, dpi = 300
)

# =======================================================================
# STEP 7 — AFTER QC SCATTER & HISTOGRAM PLOTS
# =======================================================================

df_after <- combined_after@meta.data

p1 <- ggplot(df_after, aes(nCount_RNA, nFeature_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw()

p2 <- ggplot(df_after, aes(percent.mt, nFeature_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw()

p3 <- ggplot(df_after, aes(percent.mt, nCount_RNA, color = orig.ident)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw()

p4 <- ggplot(df_after, aes(percent.mt, fill = orig.ident)) +
  geom_histogram(bins = 40, alpha = 0.7, position = "identity") +
  theme_bw()

ggsave(file.path(output_dir_after, "afterQC_scatter_nCount_vs_nFeature.jpeg"), p1, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir_after, "afterQC_scatter_percentMT_vs_nFeature.jpeg"), p2, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir_after, "afterQC_scatter_percentMT_vs_nCount.jpeg"),  p3, width = 6, height = 5, dpi = 300)
ggsave(file.path(output_dir_after, "afterQC_hist_percentMT.jpeg"),              p4, width = 6, height = 5, dpi = 300)

cat("✓ All BEFORE and AFTER QC plots saved.\n")
