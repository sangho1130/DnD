#!/gpfs/commons/home/shyoon/miniforge3/envs/scatacseq/bin/Rscript

library(ggplot2)

if(!interactive()) pdf(NULL)
args <- commandArgs(trailingOnly=TRUE)

bedFile <- args[1]
offset <- read.delim(bedFile, header = F)

ggplot(offset, aes(V5)) + 
	geom_histogram(binwidth = 1) + 
	labs(x = "Offset distance", y = "Count") + 
	theme_bw() + 
	theme(axis.text = element_text(color = "black"), 
	      axis.ticks = element_line(color = "black"),
	      panel.grid.minor = element_blank())
ggsave(gsub(bedFile, pattern = "bed", replacement = "pdf"), width = 8, height = 6, units = 'cm')


