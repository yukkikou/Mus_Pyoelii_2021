# based on david results txt
# setwd("D:/_Plasmodium_yoelii")
# after WGCNA.R

library(devtools)
library(tidyverse)
library(patchwork)
library(stringr)
library(pheatmap)

# bar of hub gene
read.delim2("figure/WGCNA/HubparasiteGene.txt", 
            header = F, sep = " ", 
            col.names = c("num","color")) %>%
  .[1:5,] %>%
  left_join(., hiCorYoeliiModule, by = "color") %>%
  ggplot(aes(x = color, y = num, fill = as.numeric(abs(cor)))) + 
    geom_bar(stat = "identity") +
    scale_fill_gradient(low = "lightblue", high = "salmon", limits=c(0.9,1)) +
    #geom_bar(stat = "identity", fill = "grey", alpha = 0.5) + 
    theme_bw() +
    theme(panel.grid = element_blank(), 
          text = element_text(size=10),
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 8,color = "black"),
          axis.text.y = element_text(size = 8,color = "black"),
          legend.position = c(0.75,0.65)) + 
    labs(x = "module", y = "Gene number", color = '') 

ggsave(file="figure/WGCNA/hubgenenum_bar.pdf",width=210/2,height=75,units="mm")

# Enrichment results
W_magenta = read_tsv("figure/WGCNA/magenta/magenta_genelist_GO.txt")
W_yellow = read_tsv("figure/WGCNA/yellow/yellow_genelist_GO.txt")
W_ltyel = read_tsv("figure/WGCNA/lightyellow/lightyellow_genelist_GO.txt")
W_turquoise = read_tsv("figure/WGCNA/turquoise/turquoise_genelist_GO.txt")
W_midnightblue = read_tsv("figure/WGCNA/midnightblue/midnightblue_genelist_GO.txt")

###******************************Enrich Term*************************###

TermAll = rbind(W_magenta,W_yellow,W_turquoise,W_ltyel,W_midnightblue) %>%
  mutate(Module = c(rep("magenta", dim(W_magenta)[1]),
                    rep("yellow", dim(W_yellow)[1]),
                    rep("turquoise", dim(W_turquoise)[1]),
                    rep("ltyel", dim(W_ltyel)[1]),
                    rep("midnightblue", dim(W_midnightblue)[1])))

TermList = rbind(W_magenta[,1:2],W_yellow[,1:2], 
                 W_turquoise[,1:2],W_ltyel[1:2],W_midnightblue[,1:2]) %>%
  distinct(Term, .keep_all = TRUE)

TermListMatrix = TermList %>% 
  mutate(magenta = TermList$Term %in% W_magenta$Term,
         yellow = TermList$Term %in% W_yellow$Term,
         ltyel = TermList$Term %in% W_ltyel$Term,
         turquoise = TermList$Term %in% W_turquoise$Term,
         midnightblue = TermList$Term %in% W_midnightblue$Term)

TermListMatrix$sum = rowSums(TermListMatrix[,3:7])
TermListMatrix = TermListMatrix[order(TermListMatrix$sum, decreasing = T),]


###******************************FILTERING*************************###
TBplus = TermListMatrix %>%
  pivot_longer(magenta:midnightblue,
               names_to = "Module",
               values_to = "own") %>%
  filter(sum > 2, own == "TRUE") %>%
  left_join(TermAll, by = c("Term" = "Term", "Module" = "Module")) %>%
  filter(grepl("KEGG_PATH", Category.x)) %>%
  KEGGtermSplit()

TB = TermListMatrix %>%
  pivot_longer(magenta:midnightblue,
               names_to = "Module",
               values_to = "own") %>%
  filter(sum > 2, own == "TRUE") %>%
  left_join(TermAll, by = c("Term" = "Term", "Module" = "Module")) %>%
  GotermSplit() %>% 
  rbind(TBplus)

TB$Module = as.factor(TB$Module)
TB$Category.x = as.factor(TB$Category.x)

TBshort = TB

TBshort[which(TB$s_Term == "defense response to Gram-positive bacterium"),]$s_Term = "defense response to Gram(+) bacterium"
TBshort[which(TB$s_Term == "positive regulation of transcription from RNA polymerase II promoter"),]$s_Term = "RNA polyII_p positive transcription regulation"
TBshort[which(TB$s_Term == "positive regulation of I-kappaB kinase/NF-kappaB signaling"),]$s_Term = "positive regulation of NF-kappaB signaling"


###******************************LAYOUT*************************###
labels <- c(GOTERM_BP_DIRECT = "BP", GOTERM_CC_DIRECT = "CC",
            GOTERM_MF_DIRECT = "MF", KEGG_PATHWAY = "K")

ggplot(TBshort,aes(s_Term, Module)) + geom_point() + 
  geom_point(aes(size=log2(`Fold Enrichment`),color=-1*log10(FDR))) + 
  facet_grid(cols = vars(Category.x), scales = "free", space = "free",
             labeller=labeller(Category.x = labels)) +
  scale_colour_gradient(low="lightsteelblue",high="lightcoral") + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=35) )+
  labs(color=expression(-log[10](FDR)),size = expression(log[2](`Fold Enrichment`)), x="",y="") + 
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "snow2"),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(color = "black"),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 5,angle = 45, hjust = 1, vjust = 1),
        text = element_text(size=10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 8, color = "black"))

ggsave(file="figure/WGCNA/heatmap.pdf",width=220,height=75*1.2,units="mm")

###******************************OUTPUT*************************###
write.csv(TB, "supplementary/multiModule_term.csv")
write.csv(TermAll, "supplementary/AllModule_term.csv")

###******************************FUNCTION*************************###
# have been ordered
GotermSplit = function(W_res){
  BP = filter(W_res, str_detect(Category.x, "GOTERM_"))
  #BP = filter(W_res, Category == Category)
  BP[, c("s_ID", "s_Term")] = str_split_fixed(BP$Term, "~", 2)
  BP
}

KEGGtermSplit = function(W_res){
  KEGG = filter(W_res, Category.x == "KEGG_PATHWAY")
  KEGG[, c("s_ID", "s_Term")] <- str_split_fixed(KEGG$Term, ":", 2)
  KEGG
}

