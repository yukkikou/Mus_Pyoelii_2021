#setwd("D:/_Plasmodium_yoelii")
#rm(list = ls())

library(tidyverse)
# pie
library(grDevices)
# PCA
library(FactoMineR)
library(factoextra)
# correlation
library(pheatmap)
# color
library(RColorBrewer)
#display.brewer.all()

# color set
col = brewer.pal(8,"Set2")

# phenotype
sampleGroup <- read.csv("expMatrix/phenotype.csv", header = TRUE)
sampleGroup$Group <- factor(sampleGroup$Group, levels = c("Control", "Treatment"))
rownames(sampleGroup) <- sampleGroup$Sample
print(sampleGroup)

# expression count matrix
circ_info <- read.csv("expMatrix/circ/circRNA_info.csv")
circ_bsj <- read.csv("expMatrix/circ/circRNA_bsj.csv")
circ_ratio <- read.csv("expMatrix/circ/circRNA_ratio.csv")
circ_de <- read.csv("expMatrix/circ/circRNA_de.tsv")


###******************************GRAPH*************************###
# classification pie of circ
pie_labels <- names(table(circ_info$circ_type))
pie_x <- table(circ_info$circ_type)
pal <- colorRampPalette(c("steelblue", "white"))
pie(pie_x, labels=pie_x, radius = 1.0, border = NA, col = pal(4))
legend("topright", pie_labels, cex=0.7, fill = pal(4))

# fall_plot
# loose thredhold for DE-circ (P < 0.05 and FC > 2)
resOrdered <- circ_de[order(circ_de$PValue),]
resSig <- subset(resOrdered, PValue < 0.05)
resSigAll <- subset(resSig, (logFC < (-1)|logFC > 1))
resSigAll$circ_id <- (rownames(resSigAll)) #28
resSigUp <- subset(resSig,logFC > 1)
length(rownames(resSigUp)) #17
resSigDown <- subset(resSig,logFC < (-1))
length(rownames(resSigDown)) #11

Siginx <- circ_de[which(circ_de$PValue<0.05),]
resRatio <- circ_ratio[which(circ_ratio$circ_id %in% Siginx$X),]
dim(resRatio)

resRatioAll <- inner_join(x = resRatio, y = Siginx, by = c("circ_id"="X"))
resCounts <- circ_bsj[which(circ_bsj$circ_id %in% Siginx$X),] %>%
  `rownames<-`(resCounts$circ_id) %>%
  select(2:9)

pre_ranked_sig_genes <- as.data.frame(resRatioAll)
pre_ranked_sig_genes <- pre_ranked_sig_genes[order(pre_ranked_sig_genes$logFC, decreasing = T),]
trend <- sapply(pre_ranked_sig_genes$logFC, function(x){if(x>0) 'up' else 'down'})
pre_ranked_sig_genes <- data.frame(pre_ranked_sig_genes, 'trend'=trend, 'rank'=1:nrow(pre_ranked_sig_genes), stringsAsFactors = F)

set.seed(12)
to_be_pointed_out <- pre_ranked_sig_genes[c(1:3, (nrow(pre_ranked_sig_genes)-2):nrow(pre_ranked_sig_genes)),] 
ggplot(pre_ranked_sig_genes, aes(x=rank, y=logFC, color=logFC))+ 
  geom_point(size=3)+ geom_hline(yintercept = c(1,-1), linetype=2, size=0.25)+ 
  geom_hline(yintercept = c(0), linetype=1, size=0.5)+ 
  geom_vline(xintercept = 17.5, linetype=2, size=0.25)+ 
  scale_color_gradient2(low="darkblue", high="darkred", mid="white", midpoint=0)+ 
  geom_point(inherit.aes = F, data=to_be_pointed_out, aes(x=rank, y=logFC), 
             size = 4, color='#F7C242')+ 
  ggrepel::geom_text_repel(inherit.aes = F, data = to_be_pointed_out, 
                           aes(x=rank, y=logFC, label=circ_id), size=3)+
  xlab('') + ylab("log2FC")+ 
  theme_bw() +
  theme(panel.grid = element_blank(), 
        text = element_text(size=10),
        panel.background = element_rect(color = 'black', fill = 'transparent'), 
        legend.key = element_rect(fill = 'transparent'),
        axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
        axis.text.y = element_text(size = 8,color = "black"),
        plot.title = element_text(hjust = 0.5),
        legend.position = c(0.9,0.68))
ggsave(file="figure/circRNA/circ_FallPlot.pdf",width=210/1.8,height=75,units="mm")

# heat map
circ_bsj %>%
  tibble::rownames_to_column(var = "circID") %>% 
  select(CASE1:CONTROL3) %>%
  as.matrix() %>% 
  t() %>% scale() %>% t() %>%
  pheatmap(color = colorRampPalette(rev(brewer.pal(n = 7, name ="RdBu")))(100),
           annotation_col = data.frame("Group" = sampleGroup$Group) %>%
             `rownames<-`(rownames(sampleGroup)), 
           annotation_legend = F,
           annotation_colors = list(
             Group = c(Control="blue3", Treatment="firebrick")))

# overlap with circAtlas

circAltas = read_tsv("../../expMatrix/circ/mm10_mm9_v2.0.txt") 
circPedia = read.csv("../../expMatrix/circ/mouse_mm10_All_circRNA.csv", header = F)


# circ gene DE
circGene_GCA39 = apply(data.frame(circ_info$gene_id), 2,
                       function(x) str_split(x, pattern = ",")) %>%
  unlist()

# circDEgene = circGene_GCA39[circGene_GCA39 %in% degSig$ensembl_gene_id]
# circDEmRNA = circGene_GCA39[circGene_GCA39 %in% degSigmRNA$ensembl_gene_id]
# circDEProPseudo = circGene_GCA39[circGene_GCA39 %in% degSigProPseudo$ensembl_gene_id]
# head(circDEgene[!(circDEgene %in% circDEmRNA)])
# # others are mainly unprocessed_pseudogene
# 
# # none
# table(circGene_GCA39 %in% degSigTEC$ensembl_gene_id)
# table(circGene_GCA39 %in% degSigLnc$ensembl_gene_id)

# circ gene classification

table(circGene_GCA39 %in% lncID)[2]
table(circGene_GCA39 %in% mRNAID)[2]
table(circGene_GCA39 %in% proPseudoID)[2]
table(circGene_GCA39 %in% unproPseudoID)[2]

length(circGene_GCA39) - (121 + 20 + 645 + 19)

#  pie
pie_labels <- c("lncRNA","mRNA","processed pseudogene","unprocessed pseudogene")
pie_x <- c(table(circGene_GCA39 %in% lncID)[2],
           table(circGene_GCA39 %in% mRNAID)[2],
           table(circGene_GCA39 %in% proPseudoID)[2],
           table(circGene_GCA39 %in% unproPseudoID)[2])
names(pie_x) = pie_labels
pal <- colorRampPalette(brewer.pal(4, "YlGnBu"))
pie(pie_x, labels=pie_x, radius = 1.0, border = NA, col = pal(4))
legend("topright", pie_labels, cex=0.7, fill = pal(4))
