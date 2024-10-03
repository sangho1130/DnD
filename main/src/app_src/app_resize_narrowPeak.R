#!/gpfs/commons/home/shyoon/miniconda3/envs/py39r43_seurat4/bin/Rscript

library(rtracklayer)

args <- commandArgs(trailingOnly=TRUE)
peak <- import(args[1], format = "narrowPeak")

peak_re <- resize(peak, width = as.numeric(args[2]), fix = "center")
peak_re <- data.frame(peak_re)

peak_w <- read.delim(args[1], header = F)
peak_w$V2 <- peak_re$start
peak_w$V3 <- peak_re$end

write.table(peak_w, paste(c(args[1], ".width", args[2]), collapse = ""), quote = F, sep = "\t", row.names = F, col.names = F)

