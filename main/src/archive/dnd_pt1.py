
# 2024-05-20; version 5.0
# 2024-09-19; finalize with SEA
# 2024-11-27; minor updates
# 2024-11-30; packaging

import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process

from main.step1 import *
from main.step2 import *
from main.step3 import *
from main.step4 import *
from main.step5 import *
from main.dnd_acc import *


def main(args):
	smt = subprocess.getoutput("which samtools")
	pcd = subprocess.getoutput("which picard")
	bcf = subprocess.getoutput("which bcftools")
	bdt = subprocess.getoutput("which bedtools")
	bdp = subprocess.getoutput("which bedops")
	mc2 = subprocess.getoutput("which macs2")
	v2b = subprocess.getoutput("which vcf2bed")

	###
	if args.Dir[-1] != "/":	args.Dir += "/"
	if not args.Output:	args.Output = args.Dir
	if args.Output[-1] != "/":	args.Output += "/"
	
	print ("### Step 1: Preprocessing ###")
	bams = flt_bam(smt, pcd, args.Dir, args.Output, args.thread, args.mapq, args.chrom, args.other, args.se, args.start, args.end)
	
	print ("### Step 2: Pileup ###")
	mpileups = mpileup_make(bcf, args.thread, args.Output, bams, args.fasta)

	print ("### Step 3: VCF & filtering ###")
	vcfs = list()
	for mpileup in mpileups:
		pu_vcf = mpileup_vcf(mpileup, args.count, args.alt)

		if not args.pass_germline:
			sort_vcf = True
			pu_vcf = flt_vcf(bdt, pu_vcf, args.gnomad, args.Output, sort_vcf)
		if args.custom:
			for customVcf in args.custom:
				sort_vcf = False
				pu_vcf = tmp_vcf(pu_vcf)
				pu_vcf = flt_vcf(bdt, pu_vcf, customVcf, args.Output, sort_vcf)

		pu_vcf = flt_vaf(pu_vcf, args.vaf) 

		vcfs.append(pu_vcf)
	print ("vcfs:", vcfs)

	print ("### Step 4: DnD edits ###")
	vcfs_flt = list()
	for vcf in vcfs:
		vcf_flt = flt_snvs(vcf, args.Output, False)
		vcfs_flt.append(vcf_flt)
	print ("vcfs_flt:", vcfs_flt)
	print ("bams:", bams)
	
	bams_flt = list()
	for bam, vcf in zip(bams, vcfs_flt):
		bam_flt = flt_alignment(bdp, smt, bam, vcf, args.Output)
		bams_flt.append(bam_flt)
	print ("bams_flt:", bams_flt)

	vcfs_dnd = list()
	for vcf in vcfs:
		vcf_ctga_flt = flt_snvs(vcf, args.Output, args.snv)
		vcfs_dnd.append(vcf_ctga_flt)
	print ("vcfs_dnd:", vcfs_dnd)

	print ("### Step 5: Peak calling ###")
	mc2dirs = list()
	rmdirs = list()

	for bam in bams:
		mc2dir = macs2_callpeak(mc2, bdt, bam, args.Output, args.gsize, args.opt, args.blacklist, args.pass_bklist, "init")
		peakFile = [mc2dir + x for x in  os.listdir(mc2dir) if x.endswith("narrowPeak")][0]
		lineCount = subprocess.getoutput("wc " + peakFile)
		if lineCount.split()[0] == 0:
			rmdirs.append(mc2dir)
			rmCmd = "rm -rf " + mc2dir
			subprocess.getoutput(rmCmd)
			continue
		mc2dirs.append(mc2dir)
	print ("mc2dirs:", mc2dirs)
	print ("removed:", rmdirs)



if __name__ == '__main__':
	parser = argparse.ArgumentParser('')
	parser.add_argument('-d', '--Dir', help = 'directory path', required = True)
	parser.add_argument('-o', '--Output', help = '[Global] (*optional) output directory path', required = False)
	parser.add_argument('--thread', help = "[Global] (*optional) number of threads; default is 4", type = str, default = "4", required = False)

	parser.add_argument('--start', type = int, help = '[Global] (*optional) the first directory index', required = False)
	parser.add_argument('--end', type = int, help = '[Global] (*optional) the last directory index', required = False)

	parser.add_argument('--mapq', help = '[Step 1] (*optional) threshold for mapping quality; default is 20', type = str, default = '20', required = False)
	parser.add_argument('--chrom', help = '[Step 1] (*optional) do not filter non-chromosome', required = False, default = False, action = 'store_true')
	parser.add_argument('--smt-other', dest = 'other', help = '[Step 1] (*optional) other filter parameters for samtools in "~"', required = False)
	parser.add_argument('--se', help = '[Step 1] (*optional) input bam is single-end', required = False, default = False, action = "store_true")

	parser.add_argument('--count', help = '[Step 2] (*optional) minimum total counts; default is 3', type = int, default = 3, required = False)
	parser.add_argument('--alt', help = '[Step 2] (*optional) minimum ALT counts; default is 2', type = int, default = 2, required = False)
	parser.add_argument('--fasta', help = '[Step 2] (*optional) genome fasta (indexed) used in alignment; e.g. cellranger/fasta/genome.fa', default = "~/DnD/genome/genome.fa", required = False)

	parser.add_argument('--gnomad', help = '[Step 3] (*optional) path to gnomAD vcf file', default = "~/DnD/filter/somatic-hg38-af-only-gnomad.hg38.chrs.vcf.gz", required = False)
	parser.add_argument('--pass-gnomad', dest = 'pass_germline', help = '[Step 3] (*optional) do not run gnomAD filering', required = False, default = False, action = 'store_true')
	parser.add_argument('--custom', help = '[Step 3] (*optional) custom vcf file(s) to filter, multiple files are accepted', nargs = "+", required = False)
	parser.add_argument('--vaf', help = '[Step 3] (*optional) filter mutations frequent than this value; e.g. 10 for 10perc; default is 10', type = int, default = 10, required = False)

	parser.add_argument('--snv', help = '[Step 4] (*optoinal) SNV patterns to search; e.g. --snv "C>T,G>A, G>C"; default is "C>T,G>A"', default = "C>T,G>A", required = False)

	parser.add_argument('--gsize', help = '[Step 5] (*optional) effective genome size for macs2 callpeak; default is hs (homo sapiens)', default = 'hs', required = False)
	parser.add_argument('--opt', help = '[Step 5] (*optional) other parameters for macs2 callpeak', type = str, required = False)
	parser.add_argument('--blacklist', help = '[Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed', required = False, default = "~/DnD/filter/hg38-blacklist.v2.bed")
	parser.add_argument('--pass-bklist', dest = 'pass_bklist', help = '[Step 5] (*optional) do not run blacklist filtering', required = False, default = False, action = 'store_true')

	args = parser.parse_args()
	main(args)


