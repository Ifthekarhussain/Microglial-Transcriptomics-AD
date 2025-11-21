# =======================================================================
# 05_sctransform_normalization.R
# Microglial Transcriptomics in Alzheimer's Disease (GSE263034)
# Author: Ifthekar Hussain
#
# Description:
#   Applies SCTransform normalization to each donor sample after
#   doublet removal. Regresses out mitochondrial percentage.
#
# Input:
#   - results/qc_filtered/seurat_list_doublet_filtered.rds
#
# Output:
#   - results/normalized/seurat_list_normalized.rds
#   - Console summary of number of samples and cells
#
# =======================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
})

set.seed(42)

# -----------------------------------------
# Define IO paths
# -----------------------------------------
input_rds   <- "results/qc_filtered/seurat_list_doublet_filtered.rds"
output_rds  <- "results/normalized/seurat_list_normalized.rds"

dir.create("results/normalized/", recursive = TRUE, showWarnings = FALSE)

# -----------------------------------------
# Load doublet-filtered samples
# -----------------------------------------
seurat_doublet_filtered <- readRDS(input_rds)

cat("Loaded", length(seurat_doublet_filtered), "samples for normalization.\n\n")

# =======================================================================
# STEP 1 — SCTransform Normalization (Per Sample)
# =======================================================================

seurat_normalized <- list()

cat("Starting SCTransform normalization...\n\n")

for (sid in names(seurat_doublet_filtered)) {

  cat("Normalizing:", sid, "\n")
  obj <- seurat_doublet_filtered[[sid]]

  # SCTransform with mitochondrial regression
  obj <- SCTransform(
    obj,
    vars.to.regress = "percent.mt",
    verbose = FALSE
  )

  seurat_normalized[[sid]] <- obj
}

cat("\n✓ SCTransform normalization completed!\n\n")

# =======================================================================
# STEP 2 — Summary
# =======================================================================

cat("=== NORMALIZATION SUMMARY ===\n")
cat("Samples normalized:", length(seurat_normalized), "\n")
cat("Total cells:", sum(sapply(seurat_normalized, ncol)), "\n")

# Save output
saveRDS(seurat_normalized, output_rds)

cat("\n✓ Saved normalized objects to:", output_rds, "\n")
