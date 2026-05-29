# Lung Squamous Cell Carcinoma (LUSC) — Bulk RNA-seq Analysis

## Project Overview

This project presents a comprehensive bulk RNA-seq analysis of **Lung Squamous Cell Carcinoma (LUSC)**, focusing on transcriptional differences between tumors originating in the **upper** and **lower lung lobes**. The study integrates differential expression analysis, functional enrichment, immune profiling, and survival analysis to explore whether anatomical lobe of origin influences gene expression and patient prognosis.

LUSC accounts for approximately 20–30% of all lung cancer cases and is strongly associated with tobacco smoking. Since both lobes share the same cell types and physiological functions, transcriptional differences are expected to be subtle — potentially more evident at the level of pathway enrichment or tumor microenvironment composition rather than large sets of differentially expressed genes.

---

## Biological Question

> **Do LUSC tumors originating in the upper lobe show a different gene expression profile compared to those in the lower lobe, and can these differences be exploited for personalized therapeutic strategies?**

---

## Analysis Pipeline

The analysis is implemented as a single R Markdown document (`LUSC_final_project.Rmd`) and is organized into the following sections:

### 1. Data Loading and Preprocessing
- Raw count matrix, sample metadata (coldata), and gene annotation (rowdata) are loaded from CSV files
- A `SummarizedExperiment` object is created and filtered to retain only tumor samples from upper and lower lobes (excluding normal tissue and ambiguous anatomical locations)

### 2. Differential Expression Analysis (DESeq2)
- A `DESeq2` object is constructed with lobe of origin as the design variable
- Low-count genes are removed using a group-size-aware filtering strategy
- Differential expression is performed comparing **upper vs lower lobe** (FDR < 0.05, |LFC| ≥ 1)
- LFC estimates are refined using **apeglm shrinkage** for more reliable effect size estimation
- **101 differentially expressed genes** were identified: 86 upregulated and 15 downregulated in the upper lobe

### 3. Visualization
- **P-value histogram** — to assess the distribution of test statistics
- **MA plots** — comparing raw and shrunken LFC estimates
- **Volcano plot** — highlighting significant genes with top hits labeled
- **Boxplots** — normalized count distributions for the top 10 DE genes across lobe groups

### 4. Functional Enrichment Analysis (clusterProfiler)
- Gene symbols are mapped to Entrez IDs using `org.Hs.eg.db`
- **Over-Representation Analysis (ORA)** is performed for GO Biological Processes (BP), GO Molecular Functions (MF), and KEGG pathways
- Results are simplified to reduce redundancy
- Key enriched themes include: amino acid metabolism, cytoskeletal remodelling, epithelial differentiation, and immune/stress response pathways

### 5. CTL–Surfaceome Correlation Analysis
- A **CTL (Cytotoxic T Lymphocyte) signature score** is computed from canonical markers: `GZMA`, `GZMB`, `GZMH`, `GZMK`, `CD8A`, `NKG7`, `PRF1`
- Expression profiles of all **surfaceome genes** (from the Human Cell Surface Proteome reference) are correlated with the CTL score using **Spearman correlation**
- Genes with |ρ| ≥ 0.75 and FDR < 0.05 are retained as CTL-correlated surface genes
- Top hits include T-cell identity markers (CD3D, CD8A, CD2), cytotoxic effectors (FASLG), chemokine receptors (CXCR3, CXCR6, CCR5), and co-inhibitory receptors (CD96, CD244)

### 6. Survival Analysis
- Kaplan–Meier curves are generated for each CTL-correlated surfaceome gene using overall survival data
- Optimal expression cutpoints are identified via `surv_cutpoint` (maximising the log-rank statistic)
- Log-rank p-values are corrected for multiple testing (Benjamini–Hochberg)
- **14 genes** reached significance at FDR < 0.05, all associated with **worse overall survival** in high-expression patients — a finding interpreted in the context of T-cell exhaustion and immunosuppressive tumor microenvironment

---

## Key Results

| Analysis | Result |
|---|---|
| DE genes identified | 101 (86 up, 15 down in upper lobe) |
| Top enriched pathways | Amino acid metabolism, cytoskeletal remodelling, immune response |
| CTL-correlated surfaceome genes (|ρ| ≥ 0.75, FDR < 0.05) | Multiple immune/lymphocyte markers |
| Survival-significant genes (FDR < 0.05) | 14 genes, all associated with worse OS |

---

## Repository Structure

```
BulkRNA/
├── LUSC_final_project.Rmd          # Main analysis script
├── LungSquamousCarcinoma_counts.csv    # Raw count matrix
├── LungSquamousCarcinoma_coldata.csv   # Sample metadata
├── LungSquamousCarcinoma_rowdata.csv   # Gene annotation
├── plots/                          # Generated figures
│   ├── DE_analysis/
│   ├── enrichment/
│   ├── ctl/
│   └── survival/
└── results/                        # Exported result tables
    ├── mapping/
    ├── enrichment/
    ├── ctl/
    └── survival/
```

---

## Requirements

### R Packages
```r
library(SummarizedExperiment)
library(DESeq2)
library(org.Hs.eg.db)
library(ggplot2)
library(clusterProfiler)
library(pheatmap)
library(readxl)
library(ggrepel)
library(survival)
library(survminer)
```

---

## How to Run

1. Clone this repository
2. Place all input CSV files in the working directory
3. Open `LUSC_final_project.Rmd` in RStudio
4. Click **Knit** or run chunks sequentially

---

## Authors

**Silvio Valli**
GitHub: [@SilvioValli](https://github.com/SilvioValli)

---

## Notes

This analysis was conducted as an exploratory study. The aim was not to perform an exhaustive biological investigation, but to provide an initial characterization of transcriptional differences between upper and lower lobe LUSC tumors. Results are intended as a starting point for future, more focused investigations.

## Disclaimer
The project description and README file were generated with the assistance of generative AI (Claude by Anthropic). The content has been reviewed and verified by the author against the original analysis documented in `LUSC_final_project.Rmd`.
