#!/gpfs/commons/home/shyoon/miniforge3/envs/scatacseq/bin/Rscript

library(ggplot2)

if(!interactive()) pdf(NULL)
args = commandArgs(trailingOnly=TRUE)

myTheme <- 
  theme_bw() + 
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.ticks = element_line(color = 'black'),
        panel.grid.minor = element_blank())

countFile_tf <- args[1]
counts_tf <- read.delim(countFile_tf)
counts_tf$SNVs <- paste(counts_tf$Ref, counts_tf$Alt, sep = ">")
counts_tf$Peak <- "Target"
counts_tf$SNVs_100bp_peak_ratio <- counts_tf$SNVs_100bp_peak/mean(counts_tf$SNVs_100bp_peak[c(1:5)])

countFile_ot <- args[2]
counts_ot <- read.delim(countFile_ot)
counts_ot$SNVs <- paste(counts_ot$Ref, counts_ot$Alt, sep = ">")
counts_ot$Peak <- "Background"
counts_ot$SNVs_100bp_peak_ratio <- counts_ot$SNVs_100bp_peak/mean(counts_ot$SNVs_100bp_peak[c(1:5)])

counts <- rbind(counts_tf, counts_ot)
counts$Peak <- factor(counts$Peak, level = c("Target", "Background"))

ggplot(counts, aes(SNVs, SNVs_100bp_peak)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Peak) +
  labs(x = "", y = "SNVs per 0.1kb per peak") +
  myTheme
ggsave(paste(c(args[3], "snvs_countNormalized.pdf"), collapse = "/"), width = 10, height = 6, units = 'cm')

ggplot(counts, aes(SNVs, SNVs_100bp_peak_ratio)) +
  geom_bar(stat = "identity") +
  facet_wrap(~Peak) +
  labs(x = "", y = "Signal to noise ratio") +
  myTheme
ggsave(paste(c(args[3], "snvs_countNormalized_SNratio.pdf"), collapse = "/"), width = 10, height = 6, units = 'cm')


