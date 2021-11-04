#setwd("D:/_Plasmodium_yoelii")
#rm(list = ls())

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


# color set
display.brewer.all()
col = brewer.pal(8,"Set2")

# phenotype
sampleGroup <- read.csv("expMatrix/phenotype.csv", header = TRUE)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("Control", "Treatment"))
rownames(sampleGroup) <- sampleGroup$Sample
print(sampleGroup)

# expression count matrix
readCount = read_csv("expMatrix/mRNA/gene_count_matrix.csv") %>%
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "gene_name"), 
                  sep = "\\|", remove = TRUE) %>%
  dplyr::select(ensembl_gene_id, rownames(sampleGroup)) %>% 
  as.data.frame()

rownames(readCount) <- readCount$ensembl_gene_id
readCount$ensembl_gene_id <- NULL
print(head(readCount, n = 3))

# ID switch
eIDs = read_csv("expMatrix/mRNA/gene_count_matrix.csv") %>%
  tidyr::separate(col = gene_id, into = c("ensembl_gene_id", "gene_name"), 
                  sep = "\\|", remove = TRUE) %>%
  dplyr::select(ensembl_gene_id, gene_name) %>% 
  as.data.frame()


# dds object
dds = DESeqDataSetFromMatrix(countData = readCount, 
                             colData = sampleGroup, 
                             design = ~ Group)


keep  = rowMeans(counts(dds)) > 1
dds = dds[keep, ]
dds = DESeq(dds)

# lfcShrink of log2FoldChange
deg = lfcShrink(dds, coef = "Group_Treatment_vs_Control", type = "apeglm") %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  dplyr::left_join(eIDs, by = "ensembl_gene_id")
degSig = dplyr::filter(deg, abs(log2FoldChange) >= 1 & padj < 0.05)

###******************************OUTPUT*************************###
write.csv(degSig, "expMatrix/mRNA/DEG_Sig.csv")

###******************************GRAPH*************************###
# normalized counts
# DESeq2 requires raw counts (without normalization)
rld <- rlog(dds)
rldt <- assay(rld) %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  dplyr::left_join(eIDs, by = "ensembl_gene_id")
print(head(rldt, n = 3))

nReadCount = counts(dds, normalized = TRUE) %>% 
  as_tibble(rownames = "ensembl_gene_id") %>% 
  dplyr::left_join(eIDs, by = "ensembl_gene_id")
print(head(nReadCount, n = 3))

# correlation
Countcor <- cor(nReadCount[,2:9])
pdf(file = "../figure/mRNA/correlation.pdf", width = 4, height = 4.5)
pheatmap(Countcor,display_numbers = TRUE,
         number_color = "black",
         fontsize_number = 8,
         color = colorRampPalette((brewer.pal(n = 9, name ="RdBu"))[c(5,9)])(200),
         border_color = NA, treeheight_row = 20, treeheight_col = 20)
dev.off()

# PCA
expr <- t(as.matrix(rldt[2:9]))
colnames(expr) = rldt$ensembl_gene_id
expr.pca <- PCA(expr, ncp = 2, scale.unit = TRUE, graph = FALSE)
#plot(expr.pca)

pca_eig1 <- round(expr.pca$eig[1,2], 2)
pca_eig2 <- round(expr.pca$eig[2,2], 2)
pca_sample <- data.frame(expr.pca$ind$coord[ ,1:2])
pca_sample <- cbind(pca_sample, sampleGroup)
pca_sample$label <- c(rep("Uninfected",3),rep("Infected",5))

ggplot(data = pca_sample, aes(x = Dim.1, y = Dim.2)) +
  geom_point(aes(color = label), size = 3) +  
  scale_color_manual(values = col) +  
  theme(panel.grid = element_blank(), 
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent')) + 
  labs(x =  paste('PCA1:', pca_eig1, '%'), y = paste('PCA2:', pca_eig2, '%'), color = '') 

ggsave(file="../figure/mRNA/pca.pdf",width=210/2,height=75,units="mm")

###******************************SPLIT*************************###
lncID = read_lines("annotation/split/Mus_musculus.GRCm39.103.lncRNA.id")
mRNAID = read_lines("annotation/split/Mus_musculus.GRCm39.103.protein_coding.id")
proPseudoID = read_lines("annotation/split/Mus_musculus.GRCm39.103.processed_pseudogene.id")
unproPseudoID = read_lines("annotation/split/Mus_musculus.GRCm39.103.unprocessed_pseudogene.id")
TECID = read_lines("annotation/split/Mus_musculus.GRCm39.103.TEC.id")
tranproPseudoID = read_lines("annotation/split/Mus_musculus.GRCm39.103.translated_unprocessed_pseudogene.id")
polyproPseudoID = read_lines("annotation/split/Mus_musculus.GRCm39.103.polymorphic_pseudogene.id")


degSigLnc = degSig %>% filter(ensembl_gene_id %in% lncID)
degSigmRNA = degSig %>% filter(ensembl_gene_id %in% mRNAID)
degSigProPseudo = degSig %>% filter(ensembl_gene_id %in% proPseudoID)
degSigTEC = degSig %>% filter(ensembl_gene_id %in% TECID)


# bar plot of DE gene in biotypes
typebar = tribble(
  ~Biotype,~DE,~Total,
  "LncRNA",dim(degSigLnc)[1],length(lncID),
  "mRNA",dim(degSigmRNA)[1],length(mRNAID),
  "Pseudogene",dim(degSigProPseudo)[1],length(proPseudoID),
  "TEC",dim(degSigTEC)[1],length(TECID),
)

other = tribble(
  ~Biotype,~DE,~Total,
  "Other",(dim(degSig)[1]-sum(typebar$DE)),(dim(readCount)[1]-sum(typebar$Total)),
  "Total",dim(degSig)[1],dim(readCount)[1]
)


rbind(typebar, other) %>%
  mutate(non_DE = Total - DE,
         DE_Ratio = DE/Total,
         factor = factor(Biotype, 
                         level = c("Total","mRNA","LncRNA","Pseudogene","TEC","Other"))) %>%
  pivot_longer(cols = c(DE, non_DE),
               names_to = "Sig",
               values_to = "number") %>%
  ggplot(aes(x = factor, y = number, fill = Sig)) + 
    geom_col(position = "stack", width = 0.65) +
    scale_fill_brewer(palette = "Set1") + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          text = element_text(size=10),
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position = c(0.8,0.8)) + 
    labs(x = "", y = "Gene number", color = '') 
ggsave(file="figure/mRNA/DE_type_bar.pdf",width=210/2,height=75,units="mm")



# volcano plot
degmRNA = deg %>% filter(ensembl_gene_id %in% mRNAID)
m1 <- volfun(degmRNA)
ggsave(file="figure/mRNA/mRNA_volcano.pdf",width=210/2,height=75,units="mm")

degLnc = deg %>% filter(ensembl_gene_id %in% lncID)
l1 <- volfun(degLnc)
ggsave(file="figure/mRNA/lncRNA_volcano.pdf",width=210/2,height=75,units="mm")


# heat map
# degSigLnc
readCount %>%
  tibble::rownames_to_column(var = "geneID") %>% 
  filter(geneID %in% degSigLnc$ensembl_gene_id) %>%
  select(COOY:XLYOY) %>%
  as.matrix() %>%
  t() %>% scale() %>% t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           annotation_col = data.frame("Group" = sampleGroup$Group) %>%
             `rownames<-`(rownames(sampleGroup)), 
           annotation_legend = F,
           annotation_colors = list(
             Group = c(Control="blue3", Treatment="firebrick")))

# degSigmRNA
readCount %>%
  tibble::rownames_to_column(var = "geneID") %>% 
  filter(geneID %in% degSigmRNA$ensembl_gene_id) %>%
  select(COOY:XLYOY) %>%
  as.matrix() %>%
  t() %>% scale() %>% t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           annotation_col = data.frame("Group" = sampleGroup$Group) %>%
             `rownames<-`(rownames(sampleGroup)), 
           annotation_legend = F,
           annotation_colors = list(
             Group = c(Control="blue3", Treatment="firebrick")))


# GO enrichment
write(degSigLnc$ensembl_gene_id, "expMatrix/mRNA/degSigLnc.txt")
write(degSigmRNA$ensembl_gene_id, "expMatrix/mRNA/degSigmRNA.txt")
# graph in GO.R

###******************************FUNCTION*************************###
volfun <- function(deg){
  for_volcano <- data.frame('logFC'=deg$log2FoldChange, 
                            'PValue'=deg$padj, 
                            'Trend' = rep('Non-sig', length(deg$padj)))
  up_sig_indices <- intersect(which(for_volcano$logFC > 1), which(for_volcano$PValue < 0.05))
  length(up_sig_indices)
  down_sig_indices <- intersect(which(for_volcano$logFC < -1), which(for_volcano$PValue < 0.05)) 
  length(down_sig_indices)
  for_volcano[up_sig_indices,'Trend'] <- 'Up'
  for_volcano[down_sig_indices, 'Trend'] <- 'Down'
  #table(for_volcano$trend)
  for_volcano$Trend <- as.factor(for_volcano$Trend)
  for_volcano$PValue <- -log10(for_volcano$PValue)
  p <- ggplot(for_volcano, aes(x=logFC, y=PValue, colour=Trend))+ 
    geom_point(size=I(1))+ 
    scale_color_manual(values = c('Non-sig'='#6E7783', 'Up'="darkred", 'Down'="darkblue"))+ 
    geom_vline(xintercept = c(1, -1), lty=2, size=I(0.8), colour = 'grey11')+ 
    geom_hline(yintercept = c(-log(x=0.05, base=10)), lty=2, size=I(0.8), colour = 'grey11')+ 
    theme_bw()+ 
    theme(panel.grid = element_blank(), 
          text = element_text(size=10),
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position = c(0.15,0.75),
          legend.text = element_text(size = 8,color = "black"),
          legend.title = element_text(size = 8,color = "black")) + 
    labs(x='log2(FoldChange)', y='-log10(adjusted p-value)')
  p
}