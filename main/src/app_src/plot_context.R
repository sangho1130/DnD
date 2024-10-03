#!/gpfs/commons/home/shyoon/miniforge3/envs/scatacseq/bin/Rscript

library(ggplot2)

if(!interactive()) pdf(NULL)
args <- commandArgs(trailingOnly=TRUE)

countFile <- args[1]
if (args[2] == 2) {
  wid <- 6
} else if (args[2] == 3) {
  wid <- 10
}

counts <- read.delim(countFile)
ggplot(counts, aes(Context, Count)) +
  geom_bar(stat = "identity") +
  labs(x = "") +
  theme_bw() + 
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.ticks = element_line(color = 'black'),
        panel.grid.minor = element_blank())
ggsave(gsub(countFile, pattern = "txt", replacement = "pdf"), width = wid, height = 6, units = 'cm')


