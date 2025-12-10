#!/usr/bin/env python3

# 

import argparse
import os
import sys
import subprocess
import importlib.resources

from main.step5 import *


def build_parser():
	parser = argparse.ArgumentParser('')
	parser.add_argument('-d', '--Dir', help = 'directory path', required = True)
	parser.add_argument('-o', '--Output', help = '[Global] (*optional) output directory path', required = False)

	parser.add_argument('--pass-peakcall', dest = 'passpc', help = '[Global] (*optional) pass peak calling and perform motif analysis ONLY; default is DO NOT PASS', required = False, default = False, action = 'store_true')
	parser.add_argument('--sample', help = '[Global, --pass-peakcall] (*optional) motif analysis for which sample: "joint" or specify sample name; default is "joint"', default = 'joint', required = False)

	parser.add_argument('--mode', help = '[Global] which mode: "sea", "homer2"; default is "sea"', nargs = '?', choices = ['sea', 'homer2'], default = 'sea', required = True)
	parser.add_argument('--fasta', help = '[Global, --mode:sea] genome fasta (indexed) used in alignment', required = True)

	parser.add_argument('--gsize', help = '[Step 5] (*optional) effective genome size for macs2 callpeak; default is hs', default = 'hs', required = False)
	parser.add_argument('--opt', help = '[Step 5] (*optional) other parameters for macs2 callpeak', type = str, required = False)
	parser.add_argument('--blacklist', help = '[Step 5] (*optional) blacklist file; default is human hg38-blacklist.v2.bed', required = False, default = None)
	parser.add_argument('--pass-bklist', dest = 'pass_bklist', help = '[Step 5] (*optional) do not run blacklist filtering', required = False, default = False, action = 'store_true')
	parser.add_argument('--motif', help = '[Step 5, --mode:sea] (*optional) motif reference; default is hocomoco v11 core human from the meme package', default = None)

	parser.add_argument('--homer-ref', dest = 'hm2ref', help = '[Step 5, --mode:homer2] (*optional) homer2 reference; default is "grch38_crgatac"', default = 'grch38_crgatac', required = False)

	return (parser)


def main(args=None):
	parser = build_parser()
	args = parser.parse_args(args)

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
		if args.blacklist == None:
			args.blacklist = str(importlib.resources.files("main")) + "/data/hg38-blacklist.v2.bed"

		mc2dir = macs2_callpeak(mc2, bdt, bams, args.Output, args.gsize, args.opt, args.blacklist, args.pass_bklist, "merged")
		print ("merged:", mc2dir)

		annotate_peaks(bdt, args.Dir, args.Output, mc2dir)

	elif args.passpc:
		if args.sample == "joint":
			mc2dir = args.Output + "step5_peaks/merged/"
		else:
			mc2dir = args.Output + "step5_peaks/" + args.sample + "/"

	if args.motif == None:
		args.motif = str(importlib.resources.files("main")) + "/data/HOCOMOCOv11_core_HUMAN_mono_meme_format.meme"
	
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
	sys.exit(main())



