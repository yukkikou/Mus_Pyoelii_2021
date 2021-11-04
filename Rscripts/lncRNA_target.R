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

# following mRNA_DEseq2.R
deLnc = nReadCount %>% 
  filter(ensembl_gene_id %in% degSigLnc$ensembl_gene_id) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  select(COOY:XLYOY)

backmRNA = nReadCount %>% 
  filter(!(ensembl_gene_id %in% degSigLnc$ensembl_gene_id)) %>%
  column_to_rownames(var = "ensembl_gene_id") %>%
  select(COOY:XLYOY)

res = NULL
for(i in 1:dim(deLnc)[1]){
  ctan = NULL
  for (j in 1:dim(backmRNA)[1]){
    corRes = cor.test(as.numeric(deLnc[i,]), as.numeric(backmRNA[j,]))
    cor = as.numeric(round(corRes$estimate, 3))
    pvalue = as.numeric(corRes$p.value)
    rs = data.frame("lncRNA"= rownames(deLnc[i,]),
                    "mRNA" = rownames(backmRNA[j,]),
                    "cor" = cor, "pvalue" = pvalue)
    ctan = rbind(ctan, rs)
  }
  res = rbind(res, ctan)
}

dim(res)

# filtering
resCor = res %>%
  mutate(adjP = p.adjust(res$pvalue, method = "BH", n = dim(res)[1]))

degCor = resCor %>%
  filter(abs(cor) >= 0.99, adjP < 0.01) %>%
  arrange(adjP)

write(unique(degCor$lncRNA), "expMatrix/lncRNA/degCorLnc.txt")
write(unique(degCor$mRNA), "expMatrix/lncRNA/degCormRNA.txt")

