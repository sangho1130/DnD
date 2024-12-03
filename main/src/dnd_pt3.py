#!/gpfs/commons/home/shyoon/miniforge3/envs/dndseq/bin/python

# 2024-05-20; version 5.0
# 2024-09-19; finalize with SEA
# 2024-11-30; packaging

import argparse
import os
import sys
import subprocess
import gzip
import random

from main.step6 import *


def main(args):
	bdt = subprocess.getoutput("which bedtools")
	v2b = subprocess.getoutput("which vcf2bed")
	if args.mode == "homer2":
		hm2 = subprocess.getoutput("which findMotifsGenome.pl")
	else:
		hm2 = False

	###
	if args.Dir[-1] != "/":	args.Dir += "/"
	if not args.Output:	args.Output = args.Dir
	if args.Output[-1] != "/":	args.Output += "/"
	
	if args.Dir[-1] != "/": args.Dir += "/"
	if not args.Output:     args.Output = args.Dir
	if args.Output[-1] != "/":      args.Output += "/"

	print ("### Step 6: Edited peaks ###")
	print ("run mode", args.mode, "\n")
	
	if args.sample == "all":
		samples = [x for x in os.listdir(args.Dir + "step5_peaks/") if x != "merged"]
		samples.sort()
	else:
		samples = [args.sample]

	for sample in samples:
		tf_peaks(args.mode, bdt, hm2, args.hm2ref, sample, args.motif, args.chipseq, args.size, args.Dir, args.Output)
		stats_peaks(args.mode, bdt, v2b, sample, args.size, args.rand, args.Dir, args.Output, args.variants)


if __name__ == '__main__':
	parser = argparse.ArgumentParser('')
	parser.add_argument('-d', '--Dir', help = 'directory path', required = True)
	parser.add_argument('-o', '--Output', help = '[Global] (*optional) output directory path', required = False)

	parser.add_argument('--size', help = '[Global] (*optional) test peak width size in bp; default is 200', type = str, default = "200", required = False)
	parser.add_argument('--rand', help = '[Global] (*optional) down-sample peaks to this number; default is 200', type = int, default = 200, required = False)
	parser.add_argument('--sample', help = '[Global] specify <sample name> or <"all"> for all samples in <step5>', required = True)
	parser.add_argument('--var', dest = 'variants', help = '[Global] (*optional) expected D&D variants; default is "C>T,G>A"', default = "C>T,G>A", required = False)

	parser.add_argument('--mode', help = '[Global] which mode: "sea", "homer2" or "chip"', nargs = '?', choices = ['chip', 'homer2', 'sea'], default = 'sea', required = True)

	parser.add_argument('--chipseq', help = '[Step 6, --mode:chip] chip-seq reference', required = False)
	parser.add_argument('--homer-ref', dest = 'hm2ref', help = '[Step 6, --mode:homer2] (*optional) homer2 reference; default is "grch38_crgatac"', default = 'grch38_crgatac', required = False)
	parser.add_argument('--motif', help = '[Step 6, --mode:homer2 or sea] <path to the homer2 motif file> for "homer2" or <TF name> for "sea", multiple arguments are supported', nargs = '+', required = False)

	args = parser.parse_args()
	main(args)


