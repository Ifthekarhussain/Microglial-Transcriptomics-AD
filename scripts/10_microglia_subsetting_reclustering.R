# =======================================================================
# 10_microglia_subsetting_reclustering.R
# Subset microglia clusters and recluster them
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
})

dir.create("results/microglia", recursive = TRUE, showWarnings = FALSE)

# Load integrated object with final clustering
combined_integrated <- readRDS("results/clustering_post_integration/combined_integrated_final_res1.0.rds")

# ---------------------------
# 1. Subset microglia clusters
# ---------------------------
microglia_clusters <- c(23, 21, 20)  # your chosen microglia cluster IDs
microglia_cells <- subset(combined_integrated, idents = microglia_clusters)

cat("Number of microglia cells:", ncol(microglia_cells), "\n")

# ---------------------------
# 2. Re-normalize in RNA assay + PCA
# ---------------------------
DefaultAssay(microglia_cells) <- "RNA"

microglia_cells <- NormalizeData(microglia_cells)
microglia_cells <- FindVariableFeatures(microglia_cells)
microglia_cells <- ScaleData(microglia_cells)
microglia_cells <- RunPCA(microglia_cells)

# Elbow plot for PCs
jpeg("results/microglia/microglia_elbow_plot.jpeg", width = 1800, height = 1200, res = 300)
ElbowPlot(microglia_cells, ndims = 30) +
  ggtitle("Elbow Plot for Microglia (30 PCs)") +
  theme_bw(base_size = 16)
dev.off()

# ---------------------------
# 3. Reclustering (15 PCs, res=1.0)
# ---------------------------
microglia_cells <- FindNeighbors(microglia_cells, dims = 1:15)
microglia_cells <- FindClusters(microglia_cells, resolution = 1.0)
microglia_cells <- RunUMAP(microglia_cells, dims = 1:15)

p_micro_umap <- DimPlot(microglia_cells, reduction = "umap", label = TRUE) +
  ggtitle("Microglia Reclustering (UMAP)")

ggsave("results/microglia/microglia_umap_clusters.jpeg", p_micro_umap,
       width = 8, height = 6, dpi = 300)

# ---------------------------
# 4. Marker detection within microglia
# ---------------------------
microglia_cells <- JoinLayers(microglia_cells)

markers <- FindAllMarkers(
  microglia_cells,
  only.pos = TRUE,
  min.pct = 0.25,
  logfc.threshold = 0.25
)

saveRDS(markers, "results/microglia/microglia_markers.rds")

top10 <- markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10, "results/microglia/microglia_top10_markers.csv", row.names = FALSE)

# ---------------------------
# 5. Rename clusters to MG states
# ---------------------------
new_ids <- c(
  "Homeostatic MG",   # 0
  "Stress MG",        # 1
  "Inflammatory MG",  # 2
  "DAM-like MG",      # 3
  "Neuron",           # 4
  "Endothelial",      # 5
  "Neuron",           # 6
  "Neuron",           # 7
  "Neuron",           # 8
  "Oligodendrocyte",  # 9
  "Fibroblast"        # 10
)

names(new_ids) <- levels(microglia_cells)
microglia_cells <- RenameIdents(microglia_cells, new_ids)

p_ann <- DimPlot(
  microglia_cells,
  reduction = "umap",
  label = TRUE,
  repel = TRUE,
  group.by = "ident"
) +
  ggtitle("Annotated Microglia & Non-Microglia Clusters")

ggsave(
  filename = "results/microglia/annotated_microglia_umap.jpeg",
  plot = p_ann,
  width = 10,
  height = 7,
  dpi = 300
)

# Save object for DEG/GO
saveRDS(microglia_cells, "results/microglia/microglia_cells_reclustered.rds")

cat("✓ Microglia subsetting & reclustering complete.\n")
