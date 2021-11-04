#setwd("D:/_Plasmodium_yoelii")
#rm(list = ls())

library(devtools)
library(tidyverse)
library(DESeq2)
library(RColorBrewer)
#PCA
library(FactoMineR)
library(factoextra)
#correlation
library(pheatmap)


# color set
col = brewer.pal(8,"Set2")

# phenotype
sampleGroup <- read.csv("expMatrix/phenotype.csv", header = TRUE)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("Control", "Treatment"))
rownames(sampleGroup) <- sampleGroup$Sample
print(sampleGroup)

# expression count matrix
hairpin = read_tsv("expMatrix/miRNA/hairpin.exp.matrix")
matrmi = read_tsv("expMatrix/miRNA/mature.exp.matrix")


# filter
HAIR_ROW_MEAN = hairpin %>%
  tibble::column_to_rownames("miR") %>%
  rowMeans() %>%
  median()

MAT_ROW_MEAN = matrmi %>%
  tibble::column_to_rownames("miR") %>%
  rowMeans() %>%
  median()

fl_hairpin = hairpin %>%
  tibble::column_to_rownames("miR") %>%
  mutate(count_mean = rowMeans(.,)) %>%
  filter(count_mean > HAIR_ROW_MEAN) %>%
  select(-count_mean)

fl_matrmi = matrmi %>%
  tibble::column_to_rownames("miR") %>%
  mutate(count_mean = rowMeans(.,)) %>%
  filter(count_mean > MAT_ROW_MEAN) %>%
  select(-count_mean)

# Different Expression
de_hairpin = dds_object(fl_hairpin)
de_matrmi = dds_object(fl_matrmi)
desSig_hairpin = de_hairpin[(de_hairpin$pvalue < 0.05 & abs(de_hairpin$log2FoldChange) > 1),]
dim(desSig_hairpin)[1] > dim(fl_hairpin)[1] * 0.05 # nothing useful

# DE-miRNA in lncRNA-seq
miID = read_lines("annotation/split/Mus_musculus.GRCm39.103.miRNA.id")
degSigmiRNA = degSig %>% filter(ensembl_gene_id %in% miID)

# verify
filter(matrmi,stringr::str_detect(miR, "142"))
filter(matrmi,stringr::str_detect(miR, "6236"))

###******************************GRAPH*************************###

# heat map
hairpin %>%
  filter(miR %in% desSig_hairpin@rownames) %>%
  tibble::column_to_rownames("miR") %>%
  select(COOY:XLYOY) %>%
  as.matrix() %>%
  t() %>% scale() %>% t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           cluster_cols = T, cluster_rows = T, annotation_names_row = T,
           treeheight_row = 10, treeheight_col = 10,
           border_color = NA, legend = F,
           annotation_col = data.frame("Group" = sampleGroup$Group) %>%
             `rownames<-`(rownames(sampleGroup)), 
           annotation_legend = F,
           annotation_colors = list(
             Group = c(Control="blue3", Treatment="firebrick")))


###******************************OUTPUT*************************###
write.csv(degSigmiRNA, "expMatrix/miRNA/DE_miRNA_Sig.csv")
write.csv(filter(de_hairpin, pvalue < 0.05), "expMatrix/miRNA/DE_hairpin_Sig.csv")


###******************************FUNCTION*************************###
dds_object <- function(readCount){
  # dds object
  dds = DESeqDataSetFromMatrix(countData = readCount, 
                               colData = sampleGroup, 
                               design = ~ Group)
  
  
  keep  = rowMeans(counts(dds)) > 1
  dds = dds[keep, ]
  dds = DESeq(dds)
  res <- results(dds)
  resOrdered <- res[order(res$padj), ]
  resOrdered
}
