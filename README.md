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

## Features
- Perform **differential expression analysis** using DESeq2  
- Map DE genes to **Homo sapiens pathways** using `org.Hs.eg.db`  
- Conduct **GSEA** using `clusterProfiler` and `msigdbr`  
- Visualize results: **NES bar plots**, **dotplots** for GO terms, **positive/negative enrichment plots**  
- Export results as **Excel files** for downstream analysis  

---

## Usage

### 1. Prepare Files
- `metadata.tsv` : Sample information (e.g., type, condition)  
- `counts.tsv` : Gene counts table from featureCounts  
### 2. Run the R Script
```r
source("GSEA.R")
```
### 3. Results
Outputs are saved in the GSEA/ folder

**Includes:**
-GSEA File.xlsx (full GSEA results)
-Dotplots for GO terms (CC, BP, MF)
-NES plots and example positive/negative enrichment plots

## License
This project is open-source under the MIT License.
You are free to use and modify it for research or educational purposes.
