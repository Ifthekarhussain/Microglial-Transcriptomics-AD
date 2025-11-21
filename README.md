# Microglial-Transcriptomics-AD

**Single-cell RNA-seq analysis of microglial transcriptomics in Alzheimer's disease**

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/Data%20Source-GSE263034-blue)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE263034)

## Overview

This repository contains a comprehensive single-cell RNA-sequencing (scRNA-seq) analysis validating key findings from **Gatenby et al. (2025)** published in *Nature Medicine*. The analysis demonstrates that amyloid-β immunization shifts microglial phenotype from metabolic stress to a neuroprotective state in Alzheimer's disease.

### Key Findings

- **nAD (non-immunized AD) microglia** exhibit metabolic stress signature with upregulated mitochondrial genes (MT-ND1-6, MT-CYB, MT-CO2, MT-ATP6)
- **iAD (immunized AD) microglia** show restored protective phenotype with upregulated APOE, TREM2, CST3, and synaptic plasticity genes
- **GO enrichment analysis** confirms metabolic dysfunction → neuroprotection recovery pathway
- Results independently validate published findings using GSE263034 dataset

---

## Repository Structure

## Repository Structure

```
Microglial-Transcriptomics-AD/
│── data/                       # Raw & processed data (ignored in .gitignore)
│   ├── raw/                    # downloaded GSE files (optional)
│   └── processed/              # Seurat objects, RDS files (optional)
│
│── scripts/                    # All R scripts for analysis
│   ├── 01_qc_filtering.R
│   ├── 02_doublet_detection.R
│   ├── 03_sctransform_normalization.R
│   ├── 04_harmony_integration.R
│   ├── 05_clustering_umap.R
│   ├── 06_deg_analysis_mast.R
│   └── 07_go_enrichment.R
│
│── results/                    # Plots (UMAPs, volcano plots, GO plots)
│   ├── qc_plots/
│   ├── umap_plots/
│   ├── volcano_plots/
│   ├── go_plots/
│   └── tables/
│
│── README.md                   # Project documentation
│── LICENSE                     # MIT license
│── .gitignore                  # Prevents uploading large data files
```


---

## Data

### Source
- **Dataset:** GSE263034 (GEO Omnibus)
- **Technology:** Single-cell RNA-sequencing (10x Genomics)
- **Samples:** 25 samples
  - NND (Non-demented controls): 6 donors
  - nAD (Non-immunized AD): 6 donors
  - iAD (Immunized AD): 13 donors
- **Cells analyzed:** 106,298 singlets (post-QC, post-doublet removal)

### Download
wget https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE263034&format=file

## Methods & Workflow

### 1. Quality Control (QC) Filtering
- Filtered cells by:
  - nFeature_RNA: median ± 2 MAD per sample
  - nCount_RNA: median ± 3 MAD per sample
  - percent.mt > 20%
- **Result:** 106,298 high-quality cells retained (from 112,196)
- Cells removed: ~5,898 doublets + low-quality cells

### 2. Doublet Detection
- Used **scDblfinder** (per-sample detection, high-confidence mode)
- Removed predicted doublets while retaining singlets
- **Result:** 94.7% singlets retained (5,898 doublets removed)

### 3. Normalization
- Applied **SCTransform** normalization
- Regressed out mitochondrial percentage (percent.mt)
- Stabilized variance across genes
- Selected highly variable genes (HVGs) for downstream analysis

### 4. Integration & Batch Correction
- **Method:** Harmony (batch correction)
- Combined all 25 samples from different donors
- **PCA:** Computed 50 PCs, selected 30 based on elbow plot
- Corrected for donor-specific batch effects
- **Result:** Samples now cluster by biology, not batch

### 5. Clustering & Cell Type Identification
- **Algorithm:** Seurat-based clustering (Louvain method)
- **Dimensionality:** UMAP visualization
- **Clusters identified:** 27 distinct cell clusters
- Filtered to microglial populations (g23, g21, g20) for differential expression

### 6. Differential Expression Analysis
- **Method:** MAST (Model-based Analysis of Single-cell Transcriptomics)
- **Statistical test:** Two-sample comparison with covariate adjustment
- **Comparisons performed:**
  - **nAD vs NND:** Disease stress signature in non-immunized AD
  - **iAD vs nAD:** Immunization rescue effect
- **Significance thresholds:** 
  - log2 Fold Change > 0.5 (or < -0.5)
  - Adjusted p-value < 0.05

### 7. Gene Ontology (GO) Enrichment Analysis
- **Package:** clusterProfiler
- **Ontology:** Biological Process (GO-BP)
- **Background:** All genes expressed in dataset
- **Enrichment method:** Hypergeometric test
- **Cutoff:** padj < 0.05
- **Result:** Identified metabolic stress and neuroprotective pathways

---

## Data Analysis Workflow Summary

Raw Data (GSE263034)
↓
QC Filtering (nFeature, nCount, percent.mt)
↓
Doublet Detection (scDblfinder)
↓
Normalization (SCTransform)
↓
Integration (Harmony)
↓
Clustering (Louvain + UMAP)
↓
Cell Type Selection (Microglial subset)
↓
Differential Expression (MAST)
↓
GO Enrichment (clusterProfiler)
↓
Validation & Interpretation


## Key Results

### Volcano Plot: nAD vs NND
**UP-regulated in nAD:**
- MT-ND1, MT-ND2, MT-ND3, MT-ND4, MT-ND5, MT-ND6
- MT-CYB, MT-CO2, MT-ATP6
- Indicates metabolic stress and compromised oxidative phosphorylation

**DOWN-regulated in nAD:**
- TMSB4X (actin regulation)
- Indicates reduced cytoskeletal support

### Volcano Plot: iAD vs nAD
**UP-regulated in iAD:**
- APOE, TREM2, CST3, S100A1 (neuroprotection)
- CKB, FAM107A (metabolic support)
- Indicates protective microglial phenotype

**DOWN-regulated in iAD:**
- HSPA1A, HSPA1B, HSP90AA1 (heat shock proteins)
- HSPH1, PTGDS
- Stress response pathways turned OFF

### GO Enrichment: nAD
- Mitochondrial organization
- ATP metabolic processes
- Macroautophagy
- **Interpretation:** nAD microglia are highly stressed, working hard to manage damage

### GO Enrichment: iAD
- Axonogenesis
- Neuron projection development
- Cognition
- Synaptic plasticity
- **Interpretation:** iAD microglia recovered neuroprotective function

---

## Usage

### Prerequisites
Install required packages
install.packages(c("Seurat", "harmony", "clusterProfiler", "org.Hs.eg.db"))
BiocManager::install(c("DESeq2", "MAST", "scDblFinder"))

---

## Key Genes & Pathways

### Metabolic Stress Signature (nAD)
| Gene | log2FC | padj | Function |
|------|--------|------|----------|
| MT-ND1 | 3.2 | 1.2e-50 | Complex I, OXPHOS |
| MT-ND2 | 2.8 | 3.4e-45 | Complex I, OXPHOS |
| MT-CYB | 2.5 | 2.1e-40 | Complex III, OXPHOS |
| MT-ATP6 | 3.5 | 8.9e-65 | ATP synthase |

### Neuroprotective Signature (iAD)
| Gene | log2FC | padj | Function |
|------|--------|------|----------|
| APOE | 1.8 | 1.2e-35 | Amyloid clearance, lipid metabolism |
| TREM2 | 1.2 | 4.5e-28 | Microglial activation, phagocytosis |
| CST3 | 2.1 | 3.7e-42 | Neuroprotection, neuroinflammation suppression |
| S100A1 | 1.5 | 2.2e-30 | Calcium signaling, synaptic plasticity |

---

## Validation

This analysis independently validates findings from:
> **Gatenby et al. (2025).** "Microglial mechanisms drive amyloid-β clearance in immunized patients with Alzheimer's disease." *Nature Medicine*, 31, 1604–1616.  
> DOI: 10.1038/s41591-025-03574-1

**Reproduced findings:**
- Mitochondrial gene upregulation in nAD microglia
-  APOE/TREM2/CST3 upregulation in immunized AD microglia
-  Heat shock protein downregulation following immunization
-  GO enrichment patterns confirm metabolic stress → neuroprotection recovery

---

## Skills & Technologies Demonstrated

**Bioinformatics Tools & Languages:**
- R programming (Seurat, Harmony, clusterProfiler, MAST)
- Single-cell RNA-seq analysis & visualization
- Statistical analysis (differential expression, GO enrichment)
- HPC cluster computing (bash, job submission)
- Git & GitHub for version control

**Computational Biology Expertise:**
- Quality control & preprocessing (QC filtering, doublet detection, normalization)
- Batch effect correction (Harmony integration)
- Cell clustering & annotation (UMAP, Louvain clustering)
- Differential expression analysis (MAST statistical testing)
- Functional annotation (Gene Ontology enrichment)

**Research Skills:**
- Literature review & paper validation
- Independent analysis & reproducibility
- Scientific interpretation & communication
- Presentation design & data visualization

---




