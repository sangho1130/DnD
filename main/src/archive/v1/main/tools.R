
library(ggplot2)
library(reshape2)

if(!interactive()) pdf(NULL)
args <- commandArgs(trailingOnly=TRUE)

myTheme <-
  theme_bw() +
  theme(axis.text = element_text(color = 'black'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = .5),
        axis.ticks = element_line(color = 'black'),
        panel.grid.minor = element_blank())


plot_offset <- function(bedFileArg) {
  offset <- read.delim(bedFileArg, header = F)

  ggplot(offset, aes(V5)) +
    geom_histogram(binwidth = 1) +
    labs(x = "Offset distance", y = "Count") +
    myTheme
  ggsave(gsub(bedFileArg, pattern = "bed", replacement = "pdf"), width = 8, height = 6, units = 'cm')
}


app_resize <- function(peakArg, sizeArg, typeArg) {
  require(rtracklayer)

  if (typeArg == "np") {
    peak <- import(peakArg, format = "narrowPeak")
  } else if (typeArg == "sm") {
    peak <- import(peakArg, format = "bed")
  }

  peak_re <- resize(peak, width = as.numeric(sizeArg), fix = "center")
  peak_re <- data.frame(peak_re)

  peak_w <- read.delim(peakArg, header = F)
  peak_w$V2 <- peak_re$start
  peak_w$V3 <- peak_re$end

  write.table(peak_w, paste(c(peakArg, ".width", sizeArg), collapse = ""), quote = F, sep = "\t", row.names = F, col.names = F)
}


plot_footprint <- function(dirPath) {
  setwd(dirPath)
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
    myTheme +
    theme(legend.position = "none")
  ggsave("footprint.pdf", units = "cm", width = 15, height = 12)
}


plot_context <- function(targetCountArg, bkgdCountArg, outPathArg) {
  countFile_tf <- targetCountArg
  counts_tf <- read.delim(countFile_tf)
  counts_tf$SNVs <- paste(counts_tf$Ref, counts_tf$Alt, sep = ">")
  counts_tf$Peak <- "Target"
  counts_tf$SNVs_peak_ratio <- counts_tf$SNVs_peak/mean(counts_tf$SNVs_peak[c(1:5)])

  countFile_ot <- bkgdCountArg
  counts_ot <- read.delim(countFile_ot)
  counts_ot$SNVs <- paste(counts_ot$Ref, counts_ot$Alt, sep = ">")
  counts_ot$Peak <- "Background"
  counts_ot$SNVs_peak_ratio <- counts_ot$SNVs_peak/mean(counts_ot$SNVs_peak[c(1:5)])

  counts <- rbind(counts_tf, counts_ot)
  counts$Peak <- factor(counts$Peak, level = c("Target", "Background"))

  ggplot(counts, aes(SNVs, SNVs_peak)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Peak) +
    labs(x = "", y = "Variants per Peak") +
    myTheme
  ggsave(paste(c(outPathArg, "dnd_countNormalized.pdf"), collapse = "/"), width = 10, height = 6, units = 'cm')

  ggplot(counts, aes(SNVs, SNVs_peak_ratio)) +
    geom_bar(stat = "identity") +
    facet_wrap(~Peak) +
    labs(x = "", y = "Signal to Noise Ratio") +
    myTheme
  ggsave(paste(c(outPathArg, "dnd_countNormalized_SNR.pdf"), collapse = "/"), width = 10, height = 6, units = 'cm')
}


extend_bed <- function(bedArg, sizeArg, varArg) {
  bed <- read.delim(bedArg, header = F)
  bed$V2 <- bed$V2 - as.numeric(sizeArg)
  bed$V3 <- bed$V3 + as.numeric(sizeArg)

  vars <- unlist(strsplit(varArg, split = "_"))
  for (var in vars) {
    ref <- unlist(strsplit(var, split = ""))[1]
    alt <- unlist(strsplit(var, split = ""))[2]
    bed_var <- subset(bed, V6 == ref & V7 == alt)
    write.table(bed_var, gsub(bedArg, pattern = "bed$", replacement = paste(c(var, "_+-", sizeArg, "bp.bed"), collapse = "")), quote = F, sep = "\t", row.names = F, col.names = F)
  }
}


ntcontext <- function(countArg, sizeArg) {
  if (sizeArg == 2) {
    wid <- 6
  } else if (sizeArg == 3) {
    wid <- 10
  }

  counts <- read.delim(countArg)
  
  ggplot(counts, aes(Context, Count)) +
    geom_bar(stat = "identity") +
    labs(x = "") +
    myTheme
  ggsave(gsub(countArg, pattern = "txt", replacement = "pdf"), width = wid, height = 6, units = 'cm')
}


###
if (args[1] == "plot_offset") {
  plot_offset(args[2])
} else if (args[1] == "app_resize") {
  app_resize(args[2], args[3], args[4])
} else if (args[1] == "plot_footprint") {
  plot_footprint(args[2])
} else if (args[1] == "plot_context") {
  plot_context(args[2], args[3], args[4])
} else if (args[1] == "extend_bed") {
  extend_bed(args[2], args[3], args[4])
} else if (args[1] == "ntcontext") {
  ntcontext(args[2], args[3])
}



