# =======================================================================
# 08_clustering_post_integration.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Tests multiple clustering resolutions on the Harmony-integrated dataset,
#   generates UMAPs for each resolution, selects final resolution (1.0),
#   and outputs final clusters.
#
# Input:
#   - results/integrated/combined_integrated.rds
#
# Output:
#   - UMAPs for each resolution
#   - Final clustering (res = 1.0)
#
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

set.seed(42)

# -----------------------------------------
# Paths
# -----------------------------------------
input_rds   <- "results/integrated/combined_integrated.rds"
plot_dir    <- "results/clustering_post_integration/"

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------
# Load integrated object
# -----------------------------------------
combined_integrated <- readRDS(input_rds)

cat("Loaded integrated dataset.\n")

# -----------------------------------------
# Test multiple resolutions
# -----------------------------------------
resolutions <- c(0.3, 0.5, 0.7, 1.0, 1.2, 1.5)
results <- list()

cat("\nTesting resolutions...\n\n")

for (res in resolutions) {
  cat("Resolution:", res, "\n")

  obj <- FindClusters(
    combined_integrated,
    resolution = res,
    verbose = FALSE
  )

  n_clust <- length(unique(obj$seurat_clusters))
  results[[as.character(res)]] <- obj

  cat(" → Clusters:", n_clust, "\n")
}

# -----------------------------------------
# Comparison Table
# -----------------------------------------
comparison <- data.frame(
  Resolution = resolutions,
  Clusters = sapply(results, \(x) length(unique(x$seurat_clusters)))
)

print(comparison)

write.csv(
  comparison,
  file.path(plot_dir, "resolution_comparison.csv"),
  row.names = FALSE
)

# -----------------------------------------
# Plot each resolution
# -----------------------------------------
cat("\nGenerating UMAPs for all resolutions...\n")

for (res in resolutions) {

  obj <- results[[as.character(res)]]
  n_clust <- length(unique(obj$seurat_clusters))

  p_res <- DimPlot(
    obj,
    reduction = "umap",
    label = TRUE,
    pt.size = 0.3
  ) +
    ggtitle(paste("Resolution", res, "→", n_clust, "clusters"))

  ggsave(
    file.path(plot_dir, paste0("resolution_", res, ".jpeg")),
    p_res, width = 7, height = 6, dpi = 300
  )
}

# -----------------------------------------
# Final resolution = 1.0 (paper used this)
# -----------------------------------------
cat("\nSelecting final resolution = 1.0\n")

combined_final <- FindClusters(
  combined_integrated,
  resolution = 1.0,
  verbose = FALSE
)

p_final <- DimPlot(
  combined_final,
  reduction = "umap",
  label = TRUE,
  pt.size = 0.5
) +
  ggtitle("Final Clustering (Resolution 1.0)")

ggsave(
  file.path(plot_dir, "final_clustering_res1.0.jpeg"),
  p_final, width = 7, height = 6, dpi = 300
)

saveRDS(
  combined_final,
  file.path(plot_dir, "combined_integrated_final_res1.0.rds")
)

cat("\n✓ Clustering (Post-Integration) complete.\n")
