#!/gpfs/commons/home/shyoon/miniconda3/envs/py39r43_seurat4/bin/Rscript

library(reshape2)
library(ggplot2)

if(!interactive()) pdf(NULL)
args <- commandArgs(trailingOnly=TRUE)

setwd(args[1])
files <- list.files("./", pattern = "_.*\\.txt$")
files <- sort(files)

data <- read.delim(files[1], sep = "\t")
rownames(data) <- data$Distance
colnames(data)[2:3] <- c("Target_1", "BG_1")

data_tf <- data[, c(1,2)]
data_bk <- data[, c(1,3)]

trial <- 2
for (file in files[-1]) {
  tmp <- read.delim(file, sep = "\t")
  
  data_tf <- cbind(data_tf, tmp[, 2])
  colnames(data_tf)[trial +1] <- paste("Target", trial, sep = "_")

  data_bk <- cbind(data_bk, tmp[, 3])
  colnames(data_bk)[trial +1] <- paste("Bkgd", trial, sep = "_")

  trial <- trial + 1
}

data_tf$Group <- "TF peaks"
data_tf$mean <- rowMeans(data_tf[, c(2:11)])
data_tf$std <- unlist( lapply( c(1:nrow(data_tf)), function(x) sd(data_tf[x, c(2:11)]) ) )

data_bk$Group <- "Background peaks"
data_bk$mean <- rowMeans(data_bk[, c(2:11)])
data_bk$std <- unlist( lapply( c(1:nrow(data_bk)), function(x) sd(data_bk[x, c(2:11)]) ) )

data <- rbind(data_tf[, c(1,12:14)], data_bk[, c(1,12:14)])
data$Distance <- as.factor(data$Distance)
data$Group <- factor(data$Group, levels = c("TF peaks", "Background peaks"))

ggplot(data, aes(x = Distance, y = mean, color = Group, group = Group)) + 
  geom_line(col = "black") + 
  geom_ribbon(aes(y = mean, ymin = mean - std, ymax = mean + std, fill = Group), alpha = .2) +
  facet_wrap(~Group, ncol = 1) +
  scale_color_manual(values = c("red2", "grey60")) +
  scale_fill_manual(values = c("red2", "grey60")) +
  scale_x_discrete(breaks = seq(from = -100, to = 100, by = 10)) + 
  labs(x = "Distance", y = "Eidt count") +
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(gsub(files[1], pattern = "_1.txt", replacement = "_line.pdf"), units = "cm", width = 15, height = 12)


snratio <- data.frame(Distance = data_tf$Distance, SN = data_tf[, "mean"]/data_bk[, "mean"])
snratio[is.nan(snratio$SN), "SN"] <- 0

ggplot(snratio, aes(x = Distance, y = SN)) + 
  geom_line() + 
  geom_smooth() +
  scale_x_continuous(breaks = seq(from = -100, to = 100, by = 10)) + 
  labs(x = "Distance", y = "SNR") +
  theme_bw() +
  theme(legend.position = "none",
        axis.ticks = element_line(color = 'black'),
        axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5))
ggsave(gsub(files[1], pattern = "_1.txt", replacement = "_SNR.pdf"), units = "cm", width = 10, height = 6)



