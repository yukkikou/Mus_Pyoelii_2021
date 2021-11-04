library(tidyverse)

list <- read_lines("../../srr_list.txt", lazy = F)
name <- unlist(strsplit(list, split = "/"))
dir <- paste(name[1],name[2],name[3],name[4],sep = "/")
subdir <- as.data.frame(sprintf("/%03d/", rep(0:99)))
names(subdir) <- "subdir"
subdir2 <- apply(subdir, 1, function(x) paste0(dir, x))
#head(subdir2)

filename <- apply(as.data.frame(rep(525:776)), 1, function(x) paste0("SRR12933", x))
#head(filename)
subfile <- filename[76:length(filename)]
#head(subfile)
fullname_600 <- paste0(subdir2, subfile)

forefile <- filename[1:75]
#head(forefile)
subdir3 <- subdir2[26:100]
fullname_525 <- paste0(subdir3, forefile)
#head(fullname_525)

fullname <- c(fullname_525, fullname_600)
#tail(fullname)

write.table(fullname, file = "SRR_url_list.txt", quote = F, row.names = F, col.names = F)

#######################fq
#fqdir <- "era-fasp@fasp.sra.ebi.ac.uk:/vol1/fastq/SRR129"
fqdir <- "/vol1/fastq/SRR129"

#demo <- "095/SRR12933595/SRR12933595.fastq.gz"

subdir_fq <- apply(subdir, 1, function(x) paste0(fqdir, x))

file_fq <- paste0(filename, ".fastq.gz")
file_fq_full <- paste(filename, file_fq, sep = "/")


subfile <- file_fq_full[76:length(file_fq_full)]
#head(subfile)
fullname_600 <- paste0(subdir_fq, subfile)

forefile <- file_fq_full[1:75]
#head(forefile)
subdir3 <- subdir_fq[26:100]
fullname_525 <- paste0(subdir3, forefile)
#head(fullname_525)

fullname <- c(fullname_525, fullname_600)

write.table(fullname, file = "fq_url_list.txt", quote = F, row.names = F, col.names = F)
