# Differential Expression Analysis: Mouse Haematopoietic Stem Cell Differentiation


## Overview
This project performs differential expression analysis on single-cell RNA-seq data from mouse haematopoietic stem cells at different stages of differentiation. The analysis compares gene expression between Haematopoietic Stem and Progenitor Cells (HSPC) and more differentiated Progenitor cells (Prog) to identify genes driving stem cell commitment. The workflow demonstrates standard bioinformatics practices including data quality control, statistical analysis using scran, and visualisation through PCA and volcano plots.


## Research Question
Which genes drive differentiation in haematopoietic stem cells as they transition from stem/progenitor states to committed progenitor cells?


## Data
- **Source**: Nestorowa et al. (2016) - single-cell RNA-seq of mouse haematopoietic stem cells.
- **Subset**: Secretome genes (genes encoding secreted proteins).
- **Sample sizes**: 
  - 701 HSPC cells (Haematopoietic Stem and Progenitor Cells).
  - 798 Prog cells (Progenitor cells).
- **Total genes analyzed**: 423 secretome genes.
- **Data format**: Log2-normalized expression values.


## Methods
### 1. Data Exploration and Quality Control
- Examined distribution of expression values across genes and cells.
- Verified data quality and sample consistency.
- Filtered for expressed genes (no genes with zero expression across all samples).

### 2. Differential Expression Analysis
- Statistical method: `scran::findMarkers()` for single-cell RNA-seq data.
- Comparison: Progenitor cells vs HSPC cells.
- Significance threshold: FDR ≤ 0.05.
- Effect size threshold: |log2 fold change| ≥ 2 for biological relevance.

### 3. Visualisation
- **Principal Component Analysis (PCA)**: Assessed separation between cell types based on gene expression profiles.
- **Volcano plots**: Visualised statistical significance vs biological effect size, highlighting key differentiation markers (e.g., Flt3).

### 4. Gene Annotation
- Retrieved gene names and descriptions from Ensembl using `biomaRt`.


## Key Results
- **282 out of 423 genes** (67%) were significantly differentially expressed between HSPC and Progenitor cells (FDR ≤ 0.05).
- Strong remodeling of the secretome during haematopoietic differentiation, indicating major changes in cell-cell signaling.
- **Top differentially expressed genes** include:
  - **Flt3** (log2FC = -3.70, FDR = 1.9×10⁻⁹⁵): Higher in HSPC, marks progenitor commitment.
  - **Ltb** (log2FC = -4.39, FDR = 1.5×10⁻¹⁴⁵): Lymphotoxin beta, involved in immune signaling.
  - **Pglyrp2** (log2FC = -4.00, FDR = 1.5×10⁻¹²⁵): Peptidoglycan recognition protein.
- PCA shows partial separation between cell types (PC1: 10.7% variance, PC2: 5.5% variance), indicating distinct but overlapping expression profiles consistent with a differentiation continuum.


## Dependencies
### R version
- R 4.3.x (or compatible version)

### Required packages
```r
# CRAN packages
install.packages(c("tidyverse", "readxl", "ggrepel"))

# Bioconductor packages
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("scran", "biomaRt"))
```

### Key packages used:
- `scran` - Single-cell differential expression analysis
- `biomaRt` - Gene annotation from Ensembl
- `tidyverse` - Data manipulation and visualization
- `ggrepel` - Text label positioning on plots


## File Structure
What's in each folder


## Usage
How to run the analysis


## References
Nestorowa, S., Hamey, F.K., Sala, B.P., Diamanti, E., Shepherd, M., Laurenti, E., Wilson, N.K., Kent, D.G. and Göttgens, B. (2016). 'A single-cell resolution map of mouse hematopoietic stem and progenitor cell differentiation'. Blood, 128(8), pp.e20–e31. doi:https://doi.org/10.1182/blood-2016-05-716480.

