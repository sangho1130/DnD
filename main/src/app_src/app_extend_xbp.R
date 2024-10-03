#!/gpfs/commons/home/shyoon/miniconda3/envs/py39r43_seurat4/bin/Rscript

args <- commandArgs(trailingOnly=TRUE)

bed <- read.delim(args[1], header = F)
bed$V2 <- bed$V2 - args[2]
bed$V3 <- bed$V3 + args[2]
bed <- subset(bed, V2 >0)

bed_ct <- subset(bed, V6 == "C" & V7 == "T")
bed_ga <- subset(bed, V6 == "G" & V7 == "A")
write.table(bed_ct, gsub(args[1], pattern = ".bed", replacement = ".CtoT_+-" + args[2] + "bp.bed"), quote = F, sep = "\t", row.names = F, col.names = F)
write.table(bed_ga, gsub(args[1], pattern = ".bed", replacement = ".GtoA_+-" + args[2] + "bp.bed"), quote = F, sep = "\t", row.names = F, col.names = F)

