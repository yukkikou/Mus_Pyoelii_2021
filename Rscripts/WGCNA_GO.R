# based on david results txt
# setwd("D:/_Plasmodium_yoelii")
# after WGCNA.R

library(devtools)
library(tidyverse)
library(patchwork)
library(stringr)
library(pheatmap)

# universal
# WGCNA/ -> combination/ 

# bar of hub gene
read.delim2("figure/combination/HubparasiteGene.txt", 
            header = F, sep = " ", 
            col.names = c("num","color")) %>%
  .[2:6,] %>%
  left_join(., hiCorYoeliiModule, by = "color") %>%
  ggplot(aes(x = reorder(color, -num), y = num, fill = as.numeric(abs(cor)))) + 
    geom_bar(stat = "identity", width = 0.5) +
    # scale_fill_gradient(low = "lightblue", high = "salmon", limits=c(0.8,1),
    #                     name = "|corrletion of modules|") +
    scale_fill_gradient(low = "lightskyblue", high = "lightgoldenrod", limits=c(0.8,1),
                      name = "|corrletion of modules|") +
    theme_bw() +
    theme(panel.grid = element_blank(), 
          text = element_text(size=8),
          panel.background = element_rect(color = 'black', fill = 'transparent'), 
          legend.key = element_rect(fill = 'transparent'),
          legend.title = element_text(size = 8, color = "black"),
          axis.text.x = element_text(hjust = 0.5, vjust = 0.5,size = 10, color = "black"),
          axis.text.y = element_text(size = 10, color = "black"),
          legend.position = c(0.8, 0.65)) + 
    labs(x = "", y = "Enriched gene number of P. yoelii", color = '')
  
ggsave(file="figure/combination/hubgenenum_bar.pdf",width=210/2,height=75,units="mm")

# Enrichment results
W_black = read_tsv("figure/combination/black/black_genelist_GO.txt")
W_brown = read_tsv("figure/combination/brown/brown_genelist_GO.txt")
W_cyan = read_tsv("figure/combination/cyan/cyan_genelist_GO.txt")
W_magenta = read_tsv("figure/combination/magenta/magenta_genelist_GO.txt")
W_pink = read_tsv("figure/combination/pink/pink_genelist_GO.txt")

###******************************Enrich Term*************************###

TermAll = rbind(W_black,W_brown,W_cyan,W_magenta,W_pink) %>%
  mutate(Module = c(rep("black", dim(W_black)[1]),
                    rep("brown", dim(W_brown)[1]),
                    rep("cyan", dim(W_cyan)[1]),
                    rep("magenta", dim(W_magenta)[1]),
                    rep("pink", dim(W_pink)[1])))

TermList = rbind(W_black[,1:2],W_brown[,1:2], 
                 W_cyan[,1:2],W_magenta[1:2],W_pink[,1:2]) %>%
  distinct(Term, .keep_all = TRUE)

TermListMatrix = TermList %>% 
  mutate(black = TermList$Term %in% W_black$Term,
         brown = TermList$Term %in% W_brown$Term,
         cyan = TermList$Term %in% W_cyan$Term,
         magenta = TermList$Term %in% W_magenta$Term,
         pink = TermList$Term %in% W_pink$Term)

TermListMatrix$sum = rowSums(TermListMatrix[,3:7])
TermListMatrix = TermListMatrix[order(TermListMatrix$sum, decreasing = T),]

###******************************FILTERING*************************###
TBplus = TermListMatrix %>%
  pivot_longer(black:pink,
               names_to = "Module",
               values_to = "own") %>%
  filter(sum > 2, own == "TRUE") %>%
  left_join(TermAll, by = c("Term" = "Term", "Module" = "Module")) %>%
  filter(grepl("KEGG_PATH", Category.x)) %>%
  KEGGtermSplit()

TB = TermListMatrix %>%
  pivot_longer(black:pink,
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
  geom_point(aes(size=log2(`Fold Enrichment`), color=-1*log10(FDR))) + 
  facet_grid(cols = vars(Category.x), scales = "free", space = "free",
             labeller=labeller(Category.x = labels)) +
  scale_colour_gradient(low="lightsteelblue",high="lightcoral") + 
  #scale_x_discrete(labels=function(x) str_wrap(x, width=35) )+
  labs(color=expression(-log[10](FDR)),size = expression(log[2](`Fold Enrichment`)), x="",y="") + 
  theme_test() +
  theme(plot.title = element_text(hjust = 0.5),
        strip.background = element_rect(fill = "snow2"),
        strip.text.x = element_text(size = 8),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 10, color = "black"),
        legend.key = element_rect(fill = 'transparent'),
        legend.key.height = unit(10, "pt"),
        legend.key.width = unit(10, "pt"),
        #axis.ticks.x = element_blank(),
        axis.text.x = element_text(size = 7, angle = 45, hjust = 1, vjust = 1),
        text = element_text(size=8, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"))

ggsave(file="figure/combination/heatmap.pdf",width=320,height=75*1.6,units="mm")

###******************************OUTPUT*************************###
write.csv(TBshort, "supplementary/multi_dual_Module_term.csv")
write.csv(TermAll, "supplementary/All_dual_Module_term.csv")

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

