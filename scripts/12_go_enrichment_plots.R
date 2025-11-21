# =======================================================================
# 12_go_enrichment_plots.R
# GO BP enrichment for DE genes and dotplots
# =======================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(stringr)
  library(dplyr)
})

dir.create("results/go", recursive = TRUE, showWarnings = FALSE)

results_nad_nnd <- readRDS("results/deg/results_nAD_vs_NND_raw.rds")
results_iad_nad <- readRDS("results/deg/results_iAD_vs_nAD_raw.rds")

# Example: take significant up-regulated genes for enrichment
sig_up_nad_nnd <- results_nad_nnd %>%
  filter(padj < 0.05, logFC > 0, !is.na(gene))

sig_up_iad_nad <- results_iad_nad %>%
  filter(padj < 0.05, logFC > 0, !is.na(gene))

gene_nad_nnd <- unique(sig_up_nad_nnd$gene)
gene_iad_nad <- unique(sig_up_iad_nad$gene)

# Map to Entrez IDs
entrez_nad_nnd <- bitr(gene_nad_nnd, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb  = org.Hs.eg.db)$ENTREZID

entrez_iad_nad <- bitr(gene_iad_nad, fromType = "SYMBOL",
                       toType = "ENTREZID",
                       OrgDb  = org.Hs.eg.db)$ENTREZID

# GO BP enrichment
go_nad_nnd <- enrichGO(
  gene         = entrez_nad_nnd,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

go_iad_nad <- enrichGO(
  gene         = entrez_iad_nad,
  OrgDb        = org.Hs.eg.db,
  ont          = "BP",
  pAdjustMethod= "BH",
  pvalueCutoff = 0.05,
  readable     = TRUE
)

# Wrap descriptions
go_nad_nnd@result$Description <- stringr::str_wrap(go_nad_nnd@result$Description, width = 28)
go_iad_nad@result$Description <- stringr::str_wrap(go_iad_nad@result$Description, width = 28)

make_go_plot <- function(go_obj, title_text) {
  dotplot(go_obj, showCategory = 15) +
    ggtitle(title_text) +
    theme_bw(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12, lineheight = 1.25),
      axis.text.x = element_text(size = 12),
      plot.title  = element_text(size = 18, face = "bold")
    )
}

p_nad_nnd <- make_go_plot(go_nad_nnd, "GO BP: nAD vs NND (Microglia)")
p_iad_nad <- make_go_plot(go_iad_nad, "GO BP: iAD vs nAD (Microglia)")

ggsave("results/go/GO_BP_nAD_vs_NND3.png", p_nad_nnd, width = 14, height = 14, dpi = 300)
ggsave("results/go/GO_BP_iAD_vs_nAD3.png", p_iad_nad, width = 14, height = 14, dpi = 300)

cat("✓ GO enrichment plots saved in results/go\n")
