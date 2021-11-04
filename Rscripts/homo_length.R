
library(tidyverse)

#args = commandArgs(T)
x = read_tsv("expMatrix/blastp/Pr_XP_022812119.1.rlt", col_names = F)
x$pr_ratio = as.numeric(apply(x, 1, function(x) min(x[3],x[6])))/as.numeric(apply(x, 1, function(x) max(x[3],x[6])))
x$ali_ratio = as.numeric(apply(x, 1, function(x) x[10])) /as.numeric(apply(x, 1, function(x) min(x[3],x[6])))







