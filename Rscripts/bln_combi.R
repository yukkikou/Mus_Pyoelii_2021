# based on blastp results txt
# setwd("D:/_Plasmodium_yoelii")

library(devtools)
library(tidyverse)
library(RColorBrewer)
#correlation
library(pheatmap)
display.brewer.all()
col = palette(brewer.pal(11, "Spectral"))


# homology gene numbers at overlap >= 90 and e < 1e20
Bn = read_table("expMatrix/blastp_combi/Final_Berghei.count", col_names = c("num", "protein"))
Cn = read_table("expMatrix/blastp_combi/Final_Chabaudi.count", col_names = c("num", "protein"))
Fn = read_table("expMatrix/blastp_combi/Final_Falciparum.count", col_names = c("num", "protein"))
Kn = read_table("expMatrix/blastp_combi/Final_Knowlesi.count", col_names = c("num", "protein"))
Mn = read_table("expMatrix/blastp_combi/Final_Malariae.count", col_names = c("num", "protein"))
On = read_table("expMatrix/blastp_combi/Final_Ovale.count", col_names = c("num", "protein"))
Vn = read_table("expMatrix/blastp_combi/Final_Vivax.count", col_names = c("num", "protein"))

j = full_join(
  full_join(
    full_join(Bn, Cn, by = "protein", suffix = c(".B", ".C")),
    full_join(Fn, Kn, by = "protein", suffix = c(".F", ".K")),
    by = "protein"
  ),
  full_join(
    full_join(Mn, On, by = "protein", suffix = c(".M", ".O")),
    Vn, by = "protein"
  ), by = "protein"
)
j[is.na(j)] = 0
j = column_to_rownames(j, var = "protein")
j$sum = apply(j, 1, sum)
j = j %>% arrange(sum)
dim(j)

# add funcion annotation
p2g = read_delim("expMatrix/blastp/map.p2g", col_names = c("protein","gene"))
highgene = read_tsv("expMatrix/blastp_combi/HighYoeliiProtein.vised.function", col_names = c("protein", "gene", "function"))

j_tmp = j %>%
  rownames_to_column(var = "protein")
j_tmp$proteinName = matrix(unlist(strsplit(j_tmp$protein,split = "\\.")),ncol=2,byrow=TRUE)[,1]

# output
proteinFun = full_join(j_tmp, highgene, by = c("proteinName" = "protein")) %>% 
  .[,2:12] %>%
  `colnames<-`(c("Bn","Cn","Fn","Kn","Mn","On","Vn","Sum","Protein","Gene","Function"))
write_tsv(proteinFun, "expMatrix/blastp_combi/Gene_homology.tsv")

## overlapping identity
# average identity
Bp = pbindn(t_file = "expMatrix/blastp_combi/Final_Berghei.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Berghei.count", species = "Berghei")
Cp = pbindn(t_file = "expMatrix/blastp_combi/Final_Chabaudi.tsv",
            c_file = "expMatrix/blastp_combi//Final_Chabaudi.count", species = "Chabaudi")
Fp = pbindn(t_file = "expMatrix/blastp_combi/Final_Falciparum.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Falciparum.count", species = "Falciparum")
Kp = pbindn(t_file = "expMatrix/blastp_combi/Final_Knowlesi.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Knowlesi.count", species = "Knowlesi")
Mp = pbindn(t_file = "expMatrix/blastp_combi/Final_Malariae.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Malariae.count", species = "Malariae")
Op = pbindn(t_file = "expMatrix/blastp_combi/Final_Ovale.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Ovale.count", species = "Ovale")
Vp = pbindn(t_file = "expMatrix/blastp_combi/Final_Vivax.tsv", 
            c_file = "expMatrix/blastp_combi/Final_Vivax.count", species = "Vivax")
# combind
pF = rbind(Bp, Cp, Fp, Kp, Mp, Op, Vp)
pF_new = pF
pF_new$protein = matrix(unlist(strsplit(pF$X1,split = "\\.")),ncol=2,byrow=TRUE)[,1]
pF_new = pF_new %>%
  left_join(highgene, by = c("protein" = "protein"))
dim(pF_new)
length(table(pF_new$X1))
length(table(pF_new$`function`))

write_tsv(pF_new, "expMatrix/blastp_combi/homology_all.tsv")

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
ggsave(file="figure/blastp_combi/highgene_YlGnBu.pdf",width=215,height=75*3.4,units="mm")

# conserved identity
Bc = cbind(file = "expMatrix/blastp_combi/Final_Berghei.tsv", species = "Berghei")
Cc = cbind(file = "expMatrix/blastp_combi/Final_Chabaudi.tsv", species = "Chabaudi")
Fc = cbind(file = "expMatrix/blastp_combi/Final_Falciparum.tsv", species = "Falciparum")
Kc = cbind(file = "expMatrix/blastp_combi/Final_Knowlesi.tsv", species = "Knowlesi")
Mc = cbind(file = "expMatrix/blastp_combi/Final_Malariae.tsv", species = "Malariae")
Oc = cbind(file = "expMatrix/blastp_combi/Final_Ovale.tsv", species = "Ovale")
Vc = cbind(file = "expMatrix/blastp_combi/Final_Vivax.tsv", species = "Vivax")

pF_single = rbind(Bc, Cc, Fc, Kc, Mc, Oc, Vc)
pF_only = pF_single
pF_only$protein = matrix(unlist(strsplit(pF_single$protein,split = "\\.")),ncol=2,byrow=TRUE)[,1]
pF_only = pF_only %>%
  left_join(highgene, by = c("protein" = "protein"))
dim(pF_only)
length(table(pF_only$protein))
length(table(pF_only$`function`))

write_tsv(pF_only, "expMatrix/blastp_combi/homology_only.tsv")

# graphic
ggplot(pF_only, aes(species, `function`))+
  #geom_tile(aes(fill="grey"),color="grey")
  geom_tile(aes(fill = pident), color="snow3")+
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
ggsave(file="figure/blastp_combi/highgeneOnly_snow3_YlGnBu.pdf",width=215,height=75*3.4,units="mm")

# conserved gene
pF_cons = pF_only %>% filter(pident > 90)
#table(rowSums(table(pF_cons$protein, pF_cons$species)) == 7)
index = (rowSums(table(pF_cons$protein, pF_cons$species)) == 7)
index = names(index[index])
pF_cons = filter(pF_cons, protein %in% index) %>%
  write_tsv(., "expMatrix/blastp_combi/conservedGene.tsv")

######################################
## function
# average
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

# only
cbind = function(file, species){
  Bc = read_delim(file, col_names = F) %>%
    .[,c(1,9:12)] %>%
    `colnames<-`(c("protein", "evalue", "length", "pident", "bitscore"))
  index = duplicated(Bc$protein)
  Bc = Bc[!index, ] %>%
    mutate(species = species)
  Bc
}

