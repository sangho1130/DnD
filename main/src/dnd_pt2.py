#!/gpfs/commons/home/shyoon/miniforge3/envs/dndseq/bin/python

# 2024-05-20; version 5.0
# 2024-09-19; finalize with SEA
# 2024-11-30; packaging

import argparse
import os
import sys
import subprocess

from main.step5 import *


def main(args):
	smt = subprocess.getoutput("which samtools")
	bdt = subprocess.getoutput("which bedtools")
	mc2 = subprocess.getoutput("which macs2")
	if args.mode == "sea":
		sea = subprocess.getoutput("which sea")
	elif args.mode == "homer2":
		hm2 = subprocess.getoutput("which findMotifsGenome.pl")

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
		mc2dir = macs2_callpeak(mc2, bdt, bams, args.Output, args.gsize, args.opt, args.blacklist, args.pass_bklist, "merged")
		print ("merged:", mc2dir)

		annotate_peaks(bdt, args.Dir, args.Output, mc2dir)

	elif args.passpc:
		mc2dir = args.Output + "step5_peaks/merged/"


	if not args.pass_bklist:
		if args.mode == "sea":
			run_memesea(bdt, sea, mc2dir, args.fasta, args.motif, "peaks_bkflt.narrowPeak", "meme_sea")
		elif args.mode == "homer2":
			run_homer2(hm2, mc2dir, args.hm2ref, "summits_bkflt.bed", "homer2")
	else:
		if args.mode == "sea":
			run_memesea(bdt, sea, mc2dir, args.fasta, args.motif, "macs2_peaks.narrowPeak", "meme_sea")
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
	parser.add_argument('--blacklist', help = '[Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed', required = False, default = "~/DnD/filter/hg38-blacklist.v2.bed")
	parser.add_argument('--pass-bklist', dest = 'pass_bklist', help = '[Step 5] (*optional) do not run blacklist filtering', required = False, default = False, action = 'store_true')
	parser.add_argument('--motif', help = '[Step 5, --mode:sea] (*optional) motif reference; default is <HOCOMOCOv11_core_HUMAN_mono_meme_format.meme>', 
				default = '/gpfs/commons/home/shyoon/Programs/meme_v5.5.5/databases/motif_databases/HUMAN/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme', required = False)
	parser.add_argument('--fasta', help = '[Step 5, --mode:sea] (*optional) genome fasta (indexed) used in alignment; e.g. cellranger/fasta/genome.fa', default = "~/DnD/genome/genome.fa", required = False)

	parser.add_argument('--homer-ref', dest = 'hm2ref', help = '[Step 5, --mode:homer2] (*optional) homer2 reference; default is "grch38_crgatac"', default = 'grch38_crgatac', required = False)

	args = parser.parse_args()
	main(args)


