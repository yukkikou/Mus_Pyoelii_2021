# based on blastp results txt
# setwd("D:/_Plasmodium_yoelii")

library(devtools)
library(tidyverse)
library(RColorBrewer)
#correlation
library(pheatmap)
display.brewer.all()
col = palette(brewer.pal(11, "Spectral"))

# e-value threshold 
fal = read_tsv("expMatrix/blastp/evalue/Final_Falciparum.txt")
ber = read_tsv("expMatrix/blastp/evalue/Final_Berghei.txt") %>%
  `colnames<-`(colnames(fal))
vix = read_tsv("expMatrix/blastp/evalue/Final_Vivax.txt") %>%
  `colnames<-`(colnames(fal))
ova = read_tsv("expMatrix/blastp/evalue/Final_Ovale.txt") %>%
  `colnames<-`(colnames(fal))
cha = read_tsv("expMatrix/blastp/evalue/Final_Chabaudi.txt") %>%
  `colnames<-`(colnames(fal))
kno = read_tsv("expMatrix/blastp/evalue/Final_Knowlesi.txt") %>%
  `colnames<-`(colnames(fal))
mal = read_tsv("expMatrix/blastp/evalue/Final_Malariae.txt") %>%
  `colnames<-`(colnames(fal))
j1 = full_join(fal[,c(1,3)], ber[,c(1,3)], by = "query_acc.ver")
j2 = full_join(vix[,c(1,3)], ova[,c(1,3)], by = "query_acc.ver")
j3 = full_join(cha[,c(1,3)], kno[,c(1,3)], by = "query_acc.ver")
j4 = full_join(j1, j2, by = "query_acc.ver")
j5 = full_join(j4, mal[,c(1,3)], by = "query_acc.ver")
j6 = full_join(j3, j5, by = "query_acc.ver") %>%
  column_to_rownames(var = "query_acc.ver") %>%
  `colnames<-`(c("fal", "ber", "vix", "ova", "cha", "kno", "mal"))
j6[is.na(j6)] = 0
htMatrix = pheatmap(j6)

# add overlap >= 90
# delete *n and used in function

# j = full_join(
#   full_join(
#     full_join(Bn, Cn, by = "gene"),
#     full_join(Fn, Kn, by = "gene"),
#     by = "gene"
#   ),
#   full_join(
#     full_join(Mn, On, by = "gene"),
#     Vn, by = "gene"
#   ), by = "gene"
# )
# j[is.na(j)] = 0
# j = column_to_rownames(j, var = "gene")
# j$sum = apply(j, 1, sum)
# j = j %>% arrange(sum)

# add funcion annotation
p2g = read_delim("expMatrix/blastp/map.p2g", col_names = c("protein","gene"))
highgene = read_tsv("expMatrix/blastp/gene.fuc", col_names = c("function", "gene", "transcript"))
pfuc = left_join(highgene, p2g, by = "gene")
write_tsv(pfuc, "expMatrix/blastp/overlap/pfuc.tsv") #remove putative

sp = matrix(unlist(strsplit(j_new$protein,split = "\\.")),ncol=2,byrow=TRUE)
j$protein = sp
proteinFun = full_join(j, pfuc, by = "protein")
write_tsv(proteinFun, "expMatrix/blastp/overlap/Gene_homology.tsv")

# identity
Bp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Berghei.tsv", c_file = "expMatrix/blastp/overlap/Final_Berghei.count", species = "Berghei")
Cp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Chabaudi.tsv", c_file = "expMatrix/blastp/overlap/Final_Chabaudi.count", species = "Chabaudi")
Fp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Falciparum.tsv", c_file = "expMatrix/blastp/overlap/Final_Falciparum.count", species = "Falciparum")
Kp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Knowlesi.tsv", c_file = "expMatrix/blastp/overlap/Final_Knowlesi.count", species = "Knowlesi")
Mp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Malariae.tsv", c_file = "expMatrix/blastp/overlap/Final_Malariae.count", species = "Malariae")
Op = pbindn(t_file = "expMatrix/blastp/overlap/Final_Ovale.tsv", c_file = "expMatrix/blastp/overlap/Final_Ovale.count", species = "Ovale")
Vp = pbindn(t_file = "expMatrix/blastp/overlap/Final_Vivax.tsv", c_file = "expMatrix/blastp/overlap/Final_Vivax.count", species = "Vivax")

pF = rbind(Bp, Cp, Fp, Kp, Mp, Op, Vp)

pF_new = pF
pfuc_new = read_tsv("expMatrix/blastp/overlap/pfuc.tsv")
sp = matrix(unlist(strsplit(pF$X1,split = "\\.")),ncol=2,byrow=TRUE)
pF_new$protein = sp[,1]
pF_new = pF_new %>%
  left_join(pfuc_new, by = c("protein" = "protein"))
dim(pF_new)
length(table(pF_new$X1))
length(table(pF_new$`function`))

# graphic
ggplot(pF_new,aes(species, `function`))+
  #geom_tile(aes(fill="grey"),color="grey")+
  geom_tile(aes(fill=identity),color="grey")+
  #scale_fill_gradient(low="lightblue",high="salmon") +
  #scale_fill_gradientn(colours = rev(brewer.pal(7, "Spectral"))) +
  scale_fill_gradientn(colours = brewer.pal(7, "YlGnBu")) +
  scale_y_discrete(labels=function(y) str_wrap(y, width=60) )+
  theme_bw() +
  labs(x="",y="") + 
  theme(panel.background = element_rect(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        panel.border = element_blank(),
        axis.text.y = element_text(color = "black"),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        axis.text.x = element_text(size = 10,angle = 45, hjust = 1, vjust = 1),
        #text = element_text(size=10, color = "black"),
        legend.text = element_text(size = 8, color = "black"),
        legend.title = element_text(size = 10, color = "black"))

ggsave(file="figure/blastp/highgene_YlGnBu.pdf",width=215,height=75*3.4,units="mm")

# function

pbindn = function(t_file, c_file, species){
  Bn = read_delim(c_file, col_names = c("n","gene")) %>% 
    mutate(species = species)
  Bp = read_tsv(t_file, col_names = F) %>%
    group_by(X1) %>%
    summarise(identity = round(mean(X11),2)) %>%
    mutate(species = species) %>%
    left_join(Bn, by = c("X1" = "gene", "species" = "species"))
  Bp
}
