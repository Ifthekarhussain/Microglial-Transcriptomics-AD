# =======================================================================
# 11_mast_deg_and_volcano.R
# Run MAST DE on microglia (nAD vs NND, iAD vs nAD) and create volcano plots
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(SingleCellExperiment)
  library(MAST)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

dir.create("results/deg", recursive = TRUE, showWarnings = FALSE)

# If you have utils_mast.R with run_mast:
# source("scripts/utils_mast.R")

microglia_cells <- readRDS("results/microglia/microglia_cells_reclustered.rds")

# ---------------------------
# 1. Subset MG only
# ---------------------------
mg_only <- subset(
  microglia_cells,
  idents = c("Homeostatic MG", "Stress MG", "Inflammatory MG", "DAM-like MG")
)

# ---------------------------
# 2. Run MAST contrasts
# ---------------------------
results_nad_nnd <- run_mast(mg_only, "NND", "nAD", "Comparison 1: nAD vs NND")
results_iad_nad <- run_mast(mg_only, "iAD", "nAD", "Comparison 2: iAD vs nAD")

# Add adjusted p-values
results_nad_nnd$padj <- p.adjust(results_nad_nnd$p, method = "BH")
results_iad_nad$padj <- p.adjust(results_iad_nad$p, method = "BH")

# Map primerid → gene
results_nad_nnd$gene <- rownames(mg_only)[match(results_nad_nnd$primerid, rownames(mg_only))]
results_iad_nad$gene <- rownames(mg_only)[match(results_iad_nad$primerid, rownames(mg_only))]

saveRDS(results_nad_nnd, "results/deg/results_nAD_vs_NND_raw.rds")
saveRDS(results_iad_nad, "results/deg/results_iAD_vs_nAD_raw.rds")

write.csv(results_nad_nnd, "results/deg/nAD_vs_NND_significant_genes.csv", row.names = FALSE)
write.csv(results_iad_nad, "results/deg/iAD_vs_nAD_significant_genes.csv", row.names = FALSE)

# ---------------------------
# 3. Helper: get top up/down genes
# ---------------------------
get_top_genes <- function(res_table, n = 20) {
  res_table <- res_table[!is.na(res_table$padj), ]
  
  up <- res_table %>%
    filter(logFC > 0) %>%
    arrange(padj) %>%
    head(n)
  
  down <- res_table %>%
    filter(logFC < 0) %>%
    arrange(padj) %>%
    head(n)
  
  list(up = up, down = down)
}

sig_nad_nnd <- results_nad_nnd %>% filter(!is.na(gene), padj < 0.05)
sig_iad_nad <- results_iad_nad %>% filter(!is.na(gene), padj < 0.05)

top20_nad_nnd <- get_top_genes(sig_nad_nnd, 20)
top20_iad_nad <- get_top_genes(sig_iad_nad, 20)

# ---------------------------
# 4. Volcano plot function (your nicer version)
# ---------------------------
plot_volcano_standard <- function(df, title, n_labels = 25,
                                  lfc_cutoff = log2(1.5),
                                  padj_cutoff = 0.05) {
  df <- df %>% filter(is.finite(logFC), is.finite(padj))
  df$logP <- -log10(df$padj)
  
  df$regulation <- case_when(
    df$padj < padj_cutoff & df$logFC >  lfc_cutoff ~ "Up-regulated",
    df$padj < padj_cutoff & df$logFC < -lfc_cutoff ~ "Down-regulated",
    TRUE ~ "Not significant"
  )
  
  df <- df %>% arrange(padj)
  df$label <- ifelse(df$regulation != "Not significant", df$primerid, "")
  
  ggplot(df, aes(x = logFC, y = logP)) +
    geom_point(aes(color = regulation), alpha = 0.8, size = 1.8) +
    scale_color_manual(values = c(
      "Up-regulated" = "red",
      "Down-regulated" = "blue",
      "Not significant" = "grey70"
    )) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff),
               linetype = "dashed", color = "black") +
    geom_hline(yintercept = -log10(padj_cutoff),
               linetype = "dashed", color = "black") +
    geom_text_repel(
      data = df %>% filter(label != "") %>% head(n_labels),
      aes(label = label),
      size = 3.3, max.overlaps = Inf
    ) +
    theme_minimal(base_size = 14) +
    labs(
      title = title,
      x = "log2 Fold Change",
      y = "-log10(FDR-adjusted p-value)",
      color = "Regulation"
    )
}

# ---------------------------
# 5. Save volcano plots
# ---------------------------
jpeg("results/deg/Volcano_nAD_vs_NND01.jpeg", width = 3000, height = 2500, res = 350)
plot_volcano_standard(results_nad_nnd, "Volcano Plot01: nAD vs NND")
dev.off()

jpeg("results/deg/Volcano_iAD_vs_nAD02.jpeg", width = 3000, height = 2500, res = 350)
plot_volcano_standard(results_iad_nad, "Volcano Plot02: iAD vs nAD")
dev.off()

cat("✓ MAST DEG + volcano plots saved in results/deg\n")
