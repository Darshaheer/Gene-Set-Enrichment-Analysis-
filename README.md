# Gene-Set-Enrichment-Analysis-
Performed Gene Set Enrichment Analysis (GSEA) on differentially expressed genes to identify enriched pathways and biological processes. Integrated DESeq2 results with clusterProfiler and msigdbr for pathway visualization in R.

## 📂 Repository Structure
Gene-Set-Enrichment-Analysis/
│── GSEA.R # R script for DESeq2 + GSEA
│── GSEA/Results/ # Results directory
│ ├── GSEA File.xlsx # GSEA results
│ ├── Ontology Plot_CC.png # Cellular Component enrichment
│ ├── Ontology Plot_BP.png # Biological Process enrichment
│ ├── Ontology Plot_MF.png # Molecular Function enrichment
│ ├── HAllmark Pathways of GSEA.png # Bar plot of NES for top pathways
│ ├── E2FPathway.png # Example positive enrichment plot
│ ├── Fatty Acid Pathway.png # Example negative enrichment plot
│── README.md # Project documentation
│── LICENSE # License file

---

This structure keeps **scripts, input files, and outputs organized**.  
All plots and results are automatically saved to the `GSEA/` folder for easy access and sharing.

---

## 🔬 Project Overview

This project combines **differential expression analysis** with **functional enrichment** to provide a complete overview of transcriptional changes in RNA-Seq experiments.

### Step 1: Differential Expression Analysis
- Input files: `metadata.tsv` (sample conditions) and `counts.tsv` (raw gene counts)  
- `DESeq2` is used to:
  - Normalize counts across samples
  - Model differential expression between experimental groups
  - Compute log2 fold changes, p-values, and adjusted p-values
- Output: Table of all genes with differential expression statistics, saved as `FINAL FILE.xlsx`

### Step 2: Gene ID Mapping
- Converts **Ensembl IDs** to gene symbols using `org.Hs.eg.db`  
- Ensures that duplicated genes are filtered by keeping the gene with the highest absolute log2 fold change  
- This step ensures proper mapping for downstream GSEA

### Step 3: Gene Set Enrichment Analysis (GSEA)
- Uses `msigdbr` to retrieve **Hallmark gene sets** for Homo sapiens  
- Creates a ranked gene list based on log2 fold change values  
- GSEA identifies pathways and gene sets enriched among upregulated or downregulated genes  
- Outputs include:
  - Enrichment scores
  - Normalized enrichment scores (NES)
  - p-values and adjusted p-values

### Step 4: GO Term Enrichment
- Uses `clusterProfiler` to perform **Gene Ontology enrichment** in three categories:
  - **Cellular Component (CC)**
  - **Biological Process (BP)**
  - **Molecular Function (MF)**
- Generates dotplots showing top enriched GO terms for each category

### Step 5: Visualization
- **Bar plots** of NES values for top pathways
- **Dotplots** for GO term enrichment
- **Positive/Negative enrichment plots** for hallmark gene sets
- All visualizations are saved as **publication-ready PNGs** in the `GSEA/` folder

---

## Features

- Handles RNA-Seq data from **any number of samples**  
- Maps **Ensembl IDs to gene symbols** for human genes  
- Computes **GSEA** and **GO enrichment** in a reproducible workflow  
- Produces **high-quality visualizations** for presentations or publications  
- Results are **organized systematically** in a dedicated `GSEA/` folder  

---

## Requirements

- R version >= 4.5.1  
- Packages required:
```r
install.packages(c("tidyverse", "dplyr", "ggplot2", "openxlsx", "magrittr"))
BiocManager::install(c(
    "DESeq2",
    "genefilter",
    "org.Hs.eg.db",
    "AnnotationDbi",
    "clusterProfiler",
    "enrichplot",
    "msigdbr"
))
```

## Usage Instructions
**1. Prepare Input Files**
- metadata.tsv should contain sample IDs and experimental conditions, e.g.,
  id	type
  sample1	normal
  sample2	treated

- counts.tsv should contain gene-level counts:
  GeneID	sample1	sample2	...
  ENSG000001234	100	200	...

**2. Run the Pipeline**
r
Copy code
source("gsea_analysis.R")

**3. Output**
  All results are saved in the GSEA/ folder, including:
  GSEA File.xlsx: Full GSEA results
  Dotplots for GO enrichment (CC, BP, MF)
  NES bar plot
  Positive and negative hallmark enrichment plots

## Example Outputs
  NES Bar Plot: Highlights pathways with the highest positive or negative enrichment
  Dotplots (CC, BP, MF): Top GO terms for enriched genes
  Positive/Negative Enrichment Plots: Visualize hallmark pathways with strongest enrichment
  Excel Files: Complete tables for downstream analysis or sharing

## Notes
  Pipeline is designed for Homo sapiens but can be adapted to other species with appropriate gene annotation
  Input counts should be gene-level and properly formatted
  Duplicate genes are automatically filtered to retain the most significant representation
  Pathways are based on Hallmark gene sets, but custom gene sets can be integrated via msigdbr
  Workflow ensures reproducibility for RNA-Seq GSEA analysis

## License
This project is open-source under the MIT License.
You are free to use, modify, and distribute it for research or educational purposes.
