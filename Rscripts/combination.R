#setwd("D:/_Plasmodium_yoelii/")

library(devtools)
library(tidyverse)
library(DESeq2)
library(ggplot2)
library(RColorBrewer)
#PCA
library(FactoMineR)
library(factoextra)
#correlation
library(pheatmap)

# dual RNA-seq matrix
CombiGeneTPM = read_csv("expMatrix/combination/combi_gene_tpm_matrix.csv")
CombiGeneCounts = read_csv("expMatrix/combination/combi_gene_count_matrix.csv")
#CombiTransTPM = read_csv("expMatrix/combination/combi_transcript_tpm_matrix.csv")
glimpse(CombiGeneTPM)
head(CombiGeneTPM)

# Mus gene
MusGene = read_csv("expMatrix/mRNA/gene_count_matrix.csv")$gene_id
head(MusGene)

# Pyoelii gene
PyoeliiGene = read_csv("expMatrix/yoelii/yoelii_gene_count_matrix.csv")$gene_id

# Mus filtering (mean(TPM) >= 1)
CombiMusGeneTPM = CombiGeneTPM %>%
  filter(gene_id %in% MusGene) %>%
  column_to_rownames(var = "gene_id") %>%
  mutate(rowmean = rowMeans(.)) %>%
  filter(rowmean >= 1) %>%
  select(-rowmean) %>%
  rownames_to_column(var = "gene_id") %>%
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "gene_name"),
                  sep = "\\|", remove = TRUE) %>%
  dplyr::select(ensembl_gene_id, rownames(sampleGroup)) %>%
  column_to_rownames(var = "ensembl_gene_id") # [1] 13316     8

# Py filtering (mean(TPM) of infected >= 1) & (control TPM == 0)
CombiPyGeneTPM = CombiGeneTPM %>%
  filter(gene_id %in% PyoeliiGene) %>%
  column_to_rownames(var = "gene_id") %>%
  filter_at(vars(starts_with("C")), all_vars((. == 0))) %>%
  mutate(rowmean = rowMeans(.)) %>%
  filter(rowmean >= 1) %>%
  select(-rowmean)  # [1] 386   8

# combinated expression matrix
head(CombiMusGeneTPM)
head(CombiPyGeneTPM)
CombiGeneFt = rbind(CombiMusGeneTPM, CombiPyGeneTPM)
dim(CombiGeneFt) # [1] 13702     8

CombiGeneCsFt = CombiGeneCounts %>%
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "gene_name"),
                  sep = "\\|", remove = TRUE) %>%
  dplyr::select(ensembl_gene_id, rownames(sampleGroup)) %>%
  filter(ensembl_gene_id %in% rownames(CombiGeneFt)) %>%
  column_to_rownames(var = "ensembl_gene_id")
  
  
