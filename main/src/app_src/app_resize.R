#!/gpfs/commons/home/shyoon/miniforge3/envs/scatacseq/bin/Rscript


resize_bed <- function(inputBed, widthSize, outputBed) {
	library(rtracklayer)

	peak <- import(inputBed, format = "BED")

	peak_re <- resize(peak, width = as.numeric(widthSize), fix = "center")
	peak_re <- data.frame(peak_re)

	peak_w <- read.delim(inputBed, header = F)
	peak_w$V2 <- peak_re$start
	peak_w$V3 <- peak_re$end

	write.table(peak_w, outputBed, quote = F, sep = "\t", row.names = F, col.names = F)

}

args <- commandArgs(trailingOnly = TRUE)
if (args[1] == '-h' | args[1] == '--help'){
    print ('Rscript app_resize.R <bed in> <width> <bed out>')
} else {
    resize_bed(args[1], args[2], args[3])
}


