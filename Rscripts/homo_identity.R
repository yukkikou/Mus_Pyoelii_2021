suppressMessages(library(tidyverse))


args = commandArgs(trailingOnly = TRUE)
file = print(args[1])

x = read_tsv(file, col_names = F, show_col_types = FALSE)

colnames(x) = c("qseqid","sseqid","qlen","qstart","qend","slen","sstart","send","evalue","length","pident","bitscore")
x$pr_ratio = round(as.numeric(apply(x, 1, function(x) min(x[3],x[6])))/as.numeric(apply(x, 1, function(x) max(x[3],x[6]))), 4)
x$ali_ratio = round(as.numeric(apply(x, 1, function(x) x[10])) /as.numeric(apply(x, 1, function(x) min(x[3],x[6]))), 4)

x_new = x %>%
filter(pr_ratio >= 0.9, ali_ratio >= 0.9)

write_tsv(x_new, paste0(x[1,1],".tsv"), col_names = F)
