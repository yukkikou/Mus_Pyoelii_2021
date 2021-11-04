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

yoeGeneCounts = read_csv("expMatrix/yoelii/yoelii_gene_count_matrix.csv")
yoeTransCounts = read_csv("expMatrix/yoelii/yoelii_transcript_count_matrix.csv")

glimpse(yoeGeneCounts)
glimpse(yoeTransCounts)

expGene = yoeGeneCounts %>%
  filter_at(vars(starts_with("C")), all_vars((. == 0))) %>%
  filter_at(vars(starts_with("X")), all_vars((. >= 5))) %>%
  column_to_rownames(var = "gene_id")

