# Loading Libraries
library(DESeq2)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(genefilter)
library(openxlsx)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(msigdbr)
library(magrittr)

# Setting Up Working Directories and Files
setwd("C:/Users/User/Downloads/")
metadata <- read.delim("metadata.tsv")
counts <- read.delim("counts.tsv")

# Differential Analysis
dseq <- DESeqDataSetFromMatrix(counts, metadata, ~type)
dseq$type <- relevel(dseq$type, "normal")
dseq <- DESeq(dseq)
res <- results(dseq)
summary(res)

resord <- as.data.frame(res)
finaltable1 <- resord[order(resord$padj), ]
write.xlsx(finaltable1, "FINAL FILE.xlsx", rowNames = TRUE)

# Preparing Directory and Files for GSEA Analysis
results_dir <- "GSEA"
if(!dir.exists(results_dir)) {
  dir.create(results_dir)
}

finaltable2 <- finaltable1
finaltable2$Gene_name <- rownames(finaltable2)
rownames(finaltable2) <- NULL

finaltable2 <- finaltable2[,c(7 , 1:6)]

dge_mapped_df <- data.frame(
  gene_syms = mapIds(
    org.Hs.eg.db,
    keys = finaltable2$Gene_name,
    column = "SYMBOL",
    keytype = "ENSEMBL",
    multiVals = "first"
  )
) %>%
  dplyr::filter(!is.na(gene_syms)) %>%
  tibble::rownames_to_column("ENSEMBL") %>%
  dplyr::inner_join(finaltable2, by = c("ENSEMBL" = "Gene_name"))

any(duplicated(dge_mapped_df$gene_syms))

dup_gene_syms <- dge_mapped_df %>%
  dplyr::filter(duplicated(gene_syms)) %>%
  dplyr::pull(gene_syms)
dge_mapped_df %>%
  dplyr::filter(gene_syms %in% dup_gene_syms) %>%
  dplyr::arrange(gene_syms)
filtered_gen_syms <- dge_mapped_df %>%
  dplyr::arrange(dplyr::desc(abs(log2FoldChange))) %>%
  dplyr::distinct(gene_syms, .keep_all = TRUE)
any(duplicated(filtered_gen_syms$gene_syms))

FC_Vector <- filtered_gen_syms$log2FoldChange
names(FC_Vector) <- filtered_gen_syms$gene_syms
FC_Vector <- sort(FC_Vector, decreasing = TRUE)
head(FC_Vector)

set.seed(2020)
msigdbr_species()
hs_hallmark <- msigdbr(
  species = "Homo sapiens",
  category = "H")
head(hs_hallmark)

# GSEA 
gsea_results <- GSEA(
  geneList = FC_Vector,
  minGSSize = 25,
  maxGSSize = 500,
  pvalueCutoff = 0.01,
  eps = 0,
  seed = TRUE,
  pAdjustMethod = "BH",
  TERM2GENE = dplyr::select(
    hs_hallmark,
    gs_name,
    gene_symbol
  )
)

head(gsea_results@result)

gsea_results_df <- data.frame(gsea_results@result)

gseaneat <- gsea_results_df %>%
  as_tibble()%>%
  arrange(desc(NES))
ggplot(gseaneat, aes(reorder(ID, NES), NES)) +
  geom_col(aes(fill=p.adjust<0.01))+
  coord_flip()+
  labs(x = "Pathways", y = "NES", title = "Hallmark Pathways of GSEA")

high_3 <- gsea_results_df %>%
  dplyr::slice_max(NES, n=3)
Positive_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_E2F_TARGETS",
  title = "HALLMARK_E2F_TARGETS",
  color.line = "blue"
)
Positive_plot

low_3 <- gsea_results_df %>%
  dplyr::slice_min(NES, n=3)
Negative_plot <- enrichplot::gseaplot(
  gsea_results,
  geneSetID = "HALLMARK_FATTY_ACID_METABOLISM",
  title = "HALLMARK_FATTY_ACID_METABOLISM",
  color.line = "red"
)
Negative_plot

write.xlsx(gsea_results_df, file.path(results_dir, "GSEA File.xlsx"), rowNames = TRUE)

ego_CC <- enrichGO(
  gene = dge_mapped_df$ENSEMBL,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "CC",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

dotplot(ego_CC, showCategory=15)

ego_BP <- enrichGO(
  gene = dge_mapped_df$ENSEMBL,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

dotplot(ego_BP, showCategory=15)

ego_MF <- enrichGO(
  gene = dge_mapped_df$ENSEMBL,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "MF",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.01,
  qvalueCutoff = 0.05
)

dotplot(ego_MF, showCategory=15)