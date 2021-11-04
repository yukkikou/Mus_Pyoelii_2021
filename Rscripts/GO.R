# based on david results txt
# setwd("D:/_Plasmodium_yoelii")
# rm(list = ls())

library(devtools)
library(tidyverse)
library(patchwork)
library(stringr)

#pvalue

# Enrichment results
enrLnc = read_tsv("expMatrix/mRNA/GoTerm/degSigLnc_GO_CC.txt") #nothing useful
enrmRNA = read_tsv("expMatrix/mRNA/GoTerm/degSigmRNA_GO_KEGG.txt")

tarmRNA = read_tsv("expMatrix/lncRNA/GoTerm/deLncTargetedmRNA.txt")

table(enrmRNA$Category)
table(tarmRNA$Category)

## enrmRNA
# have been ordered
BP = filter(enrmRNA, Category == "GOTERM_BP_DIRECT")
BP[, c("GO_ID", "GO_Term")] = str_split_fixed(BP$Term, "~", 2)  

CC = filter(enrmRNA, Category == "GOTERM_CC_DIRECT")
CC[, c("GO_ID", "GO_Term")] = str_split_fixed(CC$Term, "~", 2)  

MF = filter(enrmRNA, Category == "GOTERM_MF_DIRECT")
MF[, c("GO_ID", "GO_Term")] = str_split_fixed(MF$Term, "~", 2)  

KEGG = filter(enrmRNA, Category == "KEGG_PATHWAY")
KEGG[, c("KEGG_ID", "KEGG_Term")] <- str_split_fixed(KEGG$Term, ":", 2)

## tarmRNA
# have been ordered
tBP = filter(tarmRNA, Category == "GOTERM_BP_DIRECT")
tBP[, c("GO_ID", "GO_Term")] = str_split_fixed(tBP$Term, "~", 2)  

tCC = filter(tarmRNA, Category == "GOTERM_CC_DIRECT")
tCC[, c("GO_ID", "GO_Term")] = str_split_fixed(tCC$Term, "~", 2)  

tMF = filter(tarmRNA, Category == "GOTERM_MF_DIRECT")
tMF[, c("GO_ID", "GO_Term")] = str_split_fixed(tMF$Term, "~", 2)  

tKEGG = filter(tarmRNA, Category == "KEGG_PATHWAY")
tKEGG[, c("KEGG_ID", "KEGG_Term")] <- str_split_fixed(tKEGG$Term, ":", 2)


###******************************GRAPH*************************###
bp <- bubbleGO(BP,num=10,title="BP")
cc <-bubbleGO(CC,num=10,title="CC")
mf <-bubbleGO(MF,num=10,title="MF")
kegg <-bubbleKEGG(KEGG,num=10,title="KEGG")

bp / cc
ggsave(file="figure/mRNA/mRNA_enrichment_BP_CC.pdf",width=290/1.5,height=75*3,units="mm")

mf / kegg
ggsave(file="figure/mRNA/mRNA_enrichment_MF_KEGG.pdf",width=290/1.5,height=75*3,units="mm")

tbp <- bubbleGO(tBP,num=10,title="BP")
tcc <-bubbleGO(tCC,num=10,title="CC")
tmf <-bubbleGO(tMF,num=10,title="MF")
tkegg <-bubbleKEGG(tKEGG,num=10,title="KEGG")

tbp / tcc
ggsave(file="figure/lncRNA/tarRNA_enrichment_BP_CC.pdf",width=290/1.5,height=75*3,units="mm")

tmf / tkegg
ggsave(file="figure/lncRNA/tarRNA_enrichment_MF_KEGG.pdf",width=290/1.5,height=75*3,units="mm")


###******************************FUNCTION*************************###

bubbleGO <- function(res,num=10,title="GO term BP"){
  res <- res[c(1:num),]
  p = ggplot(res,aes(`Fold Enrichment`,GO_Term)) + geom_point() + 
    geom_point(aes(size=Count)) + 
    geom_point(aes(size=Count,color=-1*log10(FDR))) + 
    scale_y_discrete(labels=function(x) str_wrap(x, width=30) )+
    scale_colour_gradient(low="blue",high="red") + 
    labs(color=expression(-log[10](FDR)),size="Gene counts",
         x="Fold Enrichment",y="GO Term",title=title) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=12),
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 14,color = "black"),
          legend.text = element_text(size = 8,color = "black"),
          legend.title = element_text(size = 10,color = "black"))
  p
}

bubbleKEGG <- function(res,num=10,title="KEGG pathway enrichment"){
  res <- res[c(1:num),]
  p <- ggplot(res,aes(`Fold Enrichment`,KEGG_Term))+ 
    geom_point() + 
    geom_point(aes(size=Count)) +
    geom_point(aes(size=Count,color=-1*log10(FDR))) + 
    scale_colour_gradient(low="blue",high="red") +
    scale_y_discrete(labels=function(x) str_wrap(x, width=30) )+
    labs(color=expression(-log[10](FDR)),size="Gene counts",
         x="Fold Enrichment",y="KEGG Term",title=title) + 
    theme_bw() + 
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=12),
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 14,color = "black"),
          legend.text = element_text(size = 8,color = "black"),
          legend.title = element_text(size = 10,color = "black"))
  p
}
