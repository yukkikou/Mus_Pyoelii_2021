B = homo("expMatrix/blastp/overlap/Final_Berghei.tsv")
C = homo("expMatrix/blastp/overlap/Final_Chabaudi.tsv")
Fa = homo("expMatrix/blastp/overlap/Final_Falciparum.tsv")
K = homo("expMatrix/blastp/overlap/Final_Knowlesi.tsv")
M = homo("expMatrix/blastp/overlap/Final_Malariae.tsv")
O = homo("expMatrix/blastp/overlap/Final_Ovale.tsv")
V = homo("expMatrix/blastp/overlap/Final_Vivax.tsv")



Homo = c(B$X1, C$X1, Fa$X1, K$X1, M$X1, O$X1, V$X1) %>%
  table() %>% sort(., decreasing = T)

names(Homo) = matrix(unlist(strsplit(names(Homo), split = "\\.")),ncol=2,byrow=TRUE)[,1]

HighHomo = pF_new %>%
  filter(protein %in% names(Homo[Homo > 6]))


colnames(HighHomo)

ggplot(HighHomo, aes(x = `function`, y = n, fill = species)) +
  geom_bar(stat = "identity", position = "dodge") + 
  #scale_fill_manual() + 
  scale_x_discrete(labels=function(x) str_wrap(x, width=25)) +
  coord_flip() + 
  theme_classic()




# function
homo = function(file){
  B = read_tsv(file, col_names = NA) %>%
  filter(X11 > 90) %>%
  select(X1) %>%
  unique()
  B
}



