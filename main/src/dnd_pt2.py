
import argparse
import os
import sys
import subprocess


############################
### Step 5. Peak calling ###
############################

def macs2_callpeak(mcsArg, bdtArg, bamsArg, outputArg, gsizeArg, optArg, bklist, passBkArg):
	outputArg += "step5_peaks/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)
	
	sampleName = "merged"
	outputArg += sampleName + "/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	cmd = mcsArg + " callpeak"
	if not passBkArg:
		cmd += " -n __tmp_"
	else:	cmd += " -n macs2"	
	cmd += " -f BAM"
	cmd += "  -g " + gsizeArg
	cmd += " --nomodel"
	if optArg:
		cmd += " " + optArg

	cmdmc2 = cmd + " -t " + " ".join(bamsArg)
	cmdmc2 += " --outdir " + outputArg
	print (cmdmc2, "\n")
	subprocess.getoutput(cmdmc2)

	if not passBkArg:	
		cmdbt = bdtArg + " subtract -A -a " + outputArg + "__tmp__peaks.narrowPeak"
		cmdbt += " -b " + bklist
		cmdbt += " > " + outputArg + "peaks_bkflt.narrowPeak"
		print (cmdbt, "\n")
		subprocess.getoutput(cmdbt)

		openNp = open(outputArg + "peaks_bkflt.narrowPeak", "r")
		npLines = openNp.readlines()
		openNp.close()

		npLines = [x.rstrip().split() for x in npLines]
		npLines = [x[:3] + [x[3].replace("__tmp__", "")] + x[4:] for x in npLines]
		npLines = ["\t".join(x) + "\n" for x in npLines]

		openNp_w = open(outputArg + "peaks_bkflt.narrowPeak", "w")
		openNp_w.write("".join(npLines))
		openNp_w.close()

		usePeaks = [x.split()[3] for x in npLines]

		openSm = open(outputArg + "__tmp__summits.bed", "r")
		smLines = openSm.readlines()
		openSm.close()

		writeLines = list()
		for smLine in smLines:
			smLine = smLine.rstrip().split()
			smLine[3] = smLine[3].replace("__tmp__", "")
			if smLine[3] in usePeaks:
				writeLines.append("\t".join(smLine) + "\n")

		outwrite = open(outputArg + "summits_bkflt.bed", "w")
		outwrite.write("".join(writeLines))
		outwrite.close()

		cmdrm = "rm " + outputArg + "__tmp__*"
		subprocess.getoutput(cmdrm)

	return (outputArg)


def annotate_peaks(bdtArg, dirArg, outputArg, mergedOutputArg):
	outputArg += "step5_peaks/"
	
	if mergedOutputArg[-1] != "/":
		mergedOutputArg += "/"

	samples = os.listdir(dirArg + "step5_peaks/")
	samples = [x for x in samples if x != "merged"]
	samples.sort()

	mergedRes = os.listdir(mergedOutputArg)
	mergedPeaks = mergedOutputArg + [x for x in mergedRes if x.endswith("narrowPeak")][0]
	mergedSummit = mergedOutputArg + [x for x in mergedRes if x.endswith("bed")][0]

	samplePeaks = list()
	for sample in samples:
		if not os.path.exists(outputArg + sample + "_merged/"):
			os.mkdir(outputArg + sample + "_merged/")

		sampleRes = os.listdir(dirArg + "step5_peaks/" + sample)
		samplePeak = [x for x in sampleRes if x.endswith("narrowPeak")][0]
		sampleSummit = [x for x in sampleRes if x.endswith("bed")][0]

		cmd_inter = bdtArg + " intersect -wa" 
		cmd_inter += " -f 0.5 -r" # at least 50% overlap
		cmd_inter += " -a " + mergedPeaks
		cmd_inter += " -b " + dirArg + "step5_peaks/" + sample + "/" + samplePeak
		cmd_inter += " > " + outputArg + sample + "_merged/" + samplePeak + "__tmp__"
		print (cmd_inter)
		subprocess.getoutput(cmd_inter)

		## dedup
		cmd_uniq = "uniq"
		cmd_uniq += " " + outputArg + sample + "_merged/" + samplePeak + "__tmp__"
		cmd_uniq += " > " + outputArg + sample + "_merged/" + samplePeak
		print (cmd_uniq, "\n")
		subprocess.getoutput(cmd_uniq)

		subprocess.getoutput("rm " + outputArg + sample + "_merged/" + samplePeak + "__tmp__")

		if not os.path.exists(outputArg + sample + "/" + "celltype_specific"):
			os.mkdir(outputArg + sample + "/" + "celltype_specific")
		cmd_mv = "mv " + outputArg + sample + "/* " + outputArg + sample + "/celltype_specific/"
		print (cmd_mv)
		subprocess.getoutput("mv " + outputArg + sample + "/*bed " + outputArg + sample + "/celltype_specific/")
		subprocess.getoutput("mv " + outputArg + sample + "/*narrowPeak " + outputArg + sample + "/celltype_specific/")
		subprocess.getoutput("mv " + outputArg + sample + "/*fasta " + outputArg + sample + "/celltype_specific/")
		subprocess.getoutput("mv " + outputArg + sample + "/homer2 " + outputArg + sample + "/celltype_specific/")
		subprocess.getoutput("mv " + outputArg + sample + "/meme_sea " + outputArg + sample + "/celltype_specific/")	

		cmd_mv = "mv " + outputArg + sample + "_merged/* " + outputArg + sample + "/"
		print (cmd_mv)
		subprocess.getoutput(cmd_mv)
		subprocess.getoutput("rm -rf " + outputArg + sample + "_merged/")
		

		samplePeaks.append(outputArg + sample + "/" + samplePeak)
		flt_summits(outputArg + sample + "/" + samplePeak, mergedSummit)


	cmd_merge = "cat"
	cmd_merge += " " + " ".join(samplePeaks)
	cmd_merge += " > " + mergedPeaks + "__tmp__"
	print (cmd_merge, "\n")
	subprocess.getoutput(cmd_merge)
	
	if not os.path.exists(mergedOutputArg + "unfiltered"):
		os.mkdir(mergedOutputArg + "unfiltered")
	cmd_mv = "mv " + mergedOutputArg + "*narrowPeak"
	cmd_mv += " " + mergedOutputArg + "unfiltered/"
	print (cmd_mv)
	subprocess.getoutput(cmd_mv)
	cmd_mv = "mv " + mergedOutputArg + "*bed"
	cmd_mv += " " + mergedOutputArg + "unfiltered/"
	print (cmd_mv, "\n")
	subprocess.getoutput(cmd_mv)

	cmd_sort = bdtArg + " sort "
	cmd_sort += " -i " + mergedPeaks + "__tmp__"
	cmd_sort += " > " + mergedPeaks + "__dup__"
	print (cmd_sort)
	subprocess.getoutput(cmd_sort)

	cmd_rm = "rm " + mergedPeaks + "__tmp__"
	print (cmd_rm, "\n")
	subprocess.getoutput(cmd_rm)

	cmd_merge = bdtArg + " merge "
	cmd_merge += " -i " + mergedPeaks + "__dup__"
	cmd_merge += " -c " + "4,5,6,7,8,9,10"
	cmd_merge += " -o " + "first,first,first,first,first,first,first"
	cmd_merge += " > " + mergedPeaks
	print (cmd_merge)
	subprocess.getoutput(cmd_merge)

	cmd_rm = "rm " + mergedPeaks + "__dup__"
	print (cmd_rm, "\n")
	subprocess.getoutput(cmd_rm)	

	summit_merged = mergedOutputArg + "unfiltered/" + [x for x in os.listdir(mergedOutputArg + "unfiltered/") if x.endswith(".bed")][0]
	flt_summits(mergedPeaks, summit_merged)


def flt_summits(peakArg, summitArg):
	openSm = open(summitArg, "r")
	smLines = openSm.readlines()
	openSm.close()

	openPeak = open(peakArg, "r")
	peakLines = openPeak.readlines()
	openPeak.close()

	usePeaks = [x.split()[3] for x in peakLines]

	writeLines = list()
	for smLine in smLines:
		smLine = smLine.rstrip().split()
		smLine[3] = smLine[3].replace("__tmp__", "")
		if smLine[3] in usePeaks:
			writeLines.append("\t".join(smLine) + "\n")

	peakF = peakArg.split("/")[-1]
	outF = peakF.replace("peaks", "summits")
	outF = outF.replace("narrowPeak", "bed")
	outF = "/".join(peakArg.split("/")[:-1]) + "/" + outF
	outwrite = open(outF, "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


def run_homer2(hm2Arg, dirArg, hm2refArg, peakArg, outFolderArg):
	cmd_hm2 = hm2Arg
	cmd_hm2 += " " + dirArg + peakArg
	cmd_hm2 += " " + hm2refArg
	cmd_hm2 += " " + dirArg + outFolderArg
	print (cmd_hm2, "\n")
	subprocess.getoutput(cmd_hm2)


def run_memesea(bdtArg, seaArg, dirArg, motifArg, npArg, outFolderArg):
	cmd_bedtofa = bdtArg + " getfasta -nameOnly"
	cmd_bedtofa += " -fi " + "dnd_acc/genome.fa"
	cmd_bedtofa += " -fo " + dirArg + npArg.replace("narrowPeak", "fasta")
	cmd_bedtofa += " -bed " + dirArg + npArg
	print (cmd_bedtofa, "\n")
	subprocess.getoutput(cmd_bedtofa)

	cmd_sea = seaArg
	cmd_sea += " --p " + dirArg + npArg.replace("narrowPeak", "fasta")
	cmd_sea += " --m " + motifArg
	cmd_sea += " --o " + dirArg + outFolderArg
	print (cmd_sea, "\n")
	subprocess.getoutput(cmd_sea)	



def main(args):
	smt = subprocess.getoutput("which samtools")
	bdt = subprocess.getoutput("which bedtools")
	mc2 = subprocess.getoutput("which macs2")
	hm2 = subprocess.getoutput("which findMotifsGenome.pl")
	sea = subprocess.getoutput("which sea")

	###
	if args.Dir[-1] != "/":	args.Dir += "/"
	if not args.Output:	args.Output = args.Dir
	if args.Output[-1] != "/":	args.Output += "/"
	

	print ("### Step 5: Peak calling ###")
	samples = os.listdir(args.Dir + "step1_preprocess/")
	bams = list()
	for sample in samples:
		files = os.listdir(args.Dir + "step1_preprocess/" + sample)
		bamFile = args.Dir + "step1_preprocess/" + sample + "/" + [x for x in files if x.endswith(".bam")][0]
		bams.append(bamFile)
	print ("bams:", bams, "\n")

	if not args.passpc:
		mc2dir = macs2_callpeak(mc2, bdt, bams, args.Output, args.gsize, args.opt, args.blacklist, args.pass_bklist)
		print ("merged:", mc2dir)

		annotate_peaks(bdt, args.Dir, args.Output, mc2dir)

	elif args.passpc:
		mc2dir = args.Output + "step5_peaks/merged/"


	if not args.pass_bklist:
		if args.mode == "sea":
			run_memesea(bdt, sea, mc2dir, args.motif, "peaks_bkflt.narrowPeak", "meme_sea")
		elif args.mode == "homer2":
			run_homer2(hm2, mc2dir, args.hm2ref, "summits_bkflt.bed", "homer2")
	else:
		if args.mode == "sea":
			run_memesea(bdt, sea, mc2dir, args.motif, "macs2_peaks.narrowPeak", "meme_sea")
		elif args.mode == "homer2":
			run_homer2(hm2, mc2dir, args.hm2ref, "macs2_summits.bed", "homer2")

		


if __name__ == '__main__':
	parser = argparse.ArgumentParser('')
	parser.add_argument('-d', '--Dir', help = 'directory path', required = True)
	parser.add_argument('-o', '--Output', help = '[Global] (*optional) output directory path', required = False)

	parser.add_argument('--pass-peakcall', dest = 'passpc', help = '[Global] (*optional) pass peak calling and perform motif analysis ONLY; default is DO NOT PASS', required = False, default = False, action = 'store_true')

	parser.add_argument('--mode', help = '[Global] which mode: "sea", "homer2"; default is "sea"', nargs = '?', choices = ['sea', 'homer2'], default = 'sea', required = True)

	parser.add_argument('--gsize', help = '[Step 5] (*optional) effective genome size for macs2 callpeak; default is hs', default = 'hs', required = False)
	parser.add_argument('--opt', help = '[Step 5] (*optional) other parameters for macs2 callpeak', type = str, required = False)
	parser.add_argument('--blacklist', help = '[Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed', required = False, default = "dnd_acc/hg38-blacklist.v2.bed")
	parser.add_argument('--pass-bklist', dest = 'pass_bklist', help = '[Step 5] (*optional) do not run blacklist filtering', required = False, default = False, action = 'store_true')
	parser.add_argument('--motif', help = '[Step 5, --mode:sea] (*optional) motif reference; default is <HOCOMOCOv11_core_HUMAN_mono_meme_format.meme>', 
				default = 'dnd_acc/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme', required = False)
	parser.add_argument('--homer-ref', dest = 'hm2ref', help = '[Step 5, --mode:homer2] (*optional) homer2 reference; default is "grch38_crgatac"', default = 'grch38_crgatac', required = False)

	args = parser.parse_args()
	main(args)


