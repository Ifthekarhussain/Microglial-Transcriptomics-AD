# =======================================================================
# 04_doublet_detection.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Performs per-sample doublet detection using scDblFinder, removes
#   predicted doublets, summarizes results, and generates diagnostic plots.
#
# Input:
#   - seurat_list_qc.rds (from 03_qc_filtering.R)
#
# Output:
#   - seurat_list_doublet_filtered.rds
#   - Barplot of doublets removed
#   - Scatter plots after doublet removal
#   - Pie chart of overall doublet composition
#
# Saved under:
#   results/doublet_detection/
#   results/qc_plots/after_doublet/
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(dplyr)
  library(ggplot2)
})

set.seed(42)

# ---------------------------
# Define Paths
# ---------------------------
input_rds          <- "results/qc_filtered/seurat_list_qc.rds"
output_rds         <- "results/qc_filtered/seurat_list_doublet_filtered.rds"
output_dir         <- "results/doublet_detection/"
output_plot_dir    <- "results/qc_plots/after_doublet/"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(output_plot_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Load QC-filtered dataset
# ---------------------------
seurat_qc <- readRDS(input_rds)

cat("Loaded", length(seurat_qc), "samples for doublet detection.\n")

# =======================================================================
# STEP 1 — Run scDblFinder per sample
# =======================================================================

seurat_doublet_filtered <- list()

for (sid in names(seurat_qc)) {
  
  obj <- seurat_qc[[sid]]
  
  # Convert to SCE
  sce <- as.SingleCellExperiment(obj)
  
  # Run scDblFinder
  sce <- scDblFinder(sce, verbose = FALSE)
  
  # Add predictions back to Seurat
  obj$doublet_class <- sce$scDblFinder.class
  obj$doublet_score <- sce$scDblFinder.score
  
  # Keep only singlets
  obj_filtered <- subset(obj, subset = doublet_class == "singlet")
  
  seurat_doublet_filtered[[sid]] <- obj_filtered
  
  # Print summary for this sample
  cat(
    "Sample:", sid,
    "| Total:", ncol(obj),
    "| Singlets:", ncol(obj_filtered),
    "| Doublets:", ncol(obj) - ncol(obj_filtered), "\n"
  )
}

cat("\n✓ Doublet detection complete.\n")

# =======================================================================
# STEP 2 — Summary Statistics
# =======================================================================

doublet_summary <- tibble(
  Sample_ID           = names(seurat_qc),
  Cells_After_QC      = sapply(seurat_qc, ncol),
  Singlets_After_Dbl  = sapply(seurat_doublet_filtered, ncol),
  Doublets_Removed    = Cells_After_QC - Singlets_After_Dbl,
  Condition           = sapply(seurat_qc, \(x) unique(x$Condition)[1])
) %>%
  mutate(Percent_Singlet = round(100 * Singlets_After_Dbl / Cells_After_QC, 1))

cat("\n=== DOUBLET SUMMARY ===\n")
print(doublet_summary)

cat("\nOverall Singlet Retention:",
    round(100 * sum(doublet_summary$Singlets_After_Dbl) /
                 sum(doublet_summary$Cells_After_QC), 1), "%\n")

saveRDS(seurat_doublet_filtered, output_rds)

# =======================================================================
# STEP 3 — Barplot: Doublets Removed Per Sample
# =======================================================================

p_bar <- ggplot(doublet_summary, aes(x = Sample_ID, y = Doublets_Removed)) +
  geom_col(fill = "steelblue") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank()) +
  labs(
    title = "Doublets Removed Per Sample",
    x = "Samples",
    y = "Doublets Removed"
  )

ggsave(
  file.path(output_dir, "doublets_removed_barplot.jpeg"),
  p_bar, width = 8, height = 5, dpi = 300
)

# =======================================================================
# STEP 4 — Diagnostic Scatter Plots After Doublet Removal
# =======================================================================

combined_after <- Reduce(merge, seurat_doublet_filtered)
combined_after$Identity <- combined_after$orig.ident

qc_df <- combined_after@meta.data

# Scatter 1 — nCount vs nFeature
p1 <- ggplot(qc_df, aes(nCount_RNA, nFeature_RNA, color = Identity)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw() +
  labs(title = "AFTER Doublet Removal: nCount vs nFeature")

ggsave(
  file.path(output_plot_dir, "after_doublet_nCount_vs_nFeature.jpeg"),
  p1, width = 6, height = 5, dpi = 300
)

# Scatter 2 — percent.mt vs nCount
p2 <- ggplot(qc_df, aes(percent.mt, nCount_RNA, color = Identity)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw() +
  labs(title = "AFTER Doublet Removal: percent.mt vs nCount")

ggsave(
  file.path(output_plot_dir, "after_doublet_percentMT_vs_nCount.jpeg"),
  p2, width = 6, height = 5, dpi = 300
)

# Scatter 3 — percent.mt vs nFeature
p3 <- ggplot(qc_df, aes(percent.mt, nFeature_RNA, color = Identity)) +
  geom_point(alpha = 0.4, size = 0.7) +
  theme_bw() +
  labs(title = "AFTER Doublet Removal: percent.mt vs nFeature")

ggsave(
  file.path(output_plot_dir, "after_doublet_percentMT_vs_nFeature.jpeg"),
  p3, width = 6, height = 5, dpi = 300
)

# =======================================================================
# STEP 5 — Pie Chart of Overall Doublet Composition
# =======================================================================

total_singlet <- sum(doublet_summary$Singlets_After_Dbl)
total_doublet <- sum(doublet_summary$Doublets_Removed)

df_pie <- tibble(
  Class  = c("Singlets", "Doublets"),
  Count  = c(total_singlet, total_doublet)
) %>%
  mutate(
    Percent = round(Count / sum(Count) * 100, 1),
    Label = paste0(Class, " (", Percent, "%)")
  )

p_pie <- ggplot(df_pie, aes(x = "", y = Count, fill = Class)) +
  geom_col(width = 1) +
  coord_polar("y") +
  geom_label(aes(label = Label), position = position_stack(0.5), size = 6) +
  scale_fill_manual(values = c("skyblue", "coral")) +
  theme_void() +
  labs(title = "Overall Doublet Composition") +
  theme(plot.title = element_text(size = 18, hjust = 0.5))

ggsave(
  file.path(output_dir, "doublet_pie_chart.jpeg"),
  p_pie, width = 6, height = 6, dpi = 300
)

cat("\n✓ Doublet detection outputs saved successfully.\n")
