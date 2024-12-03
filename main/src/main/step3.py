
import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process
from main.dnd_acc import *


def flt_vcf(bdtArg, vcfArg, rmArg, outputArg, sortArg):
	outputArg += "step3_fltvcf/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	sampleName = vcfArg.split("/")[-2]
	outputArg += sampleName + "/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	cmd_sort = 'grep "^#" ' + vcfArg + " > " + outputArg + "__tmp__.sorted.vcf"
	cmd_sort += ' && grep -v "^#" ' + vcfArg + " | sort -V -k1,1 -k2,2n >> " + outputArg + "__tmp__.sorted.vcf"
	print (cmd_sort, "\n")
	subprocess.getoutput(cmd_sort)

	cmd_sub = bdtArg + " subtract -header"
	if sortArg:
		cmd_sub += " -sorted"
	cmd_sub += " -a " + outputArg + "__tmp__.sorted.vcf"
	cmd_sub += " -b " + rmArg
	cmd_sub += " > " + outputArg + "__tmp_subtracted__.vcf"
	print (cmd_sub, "\n")
	subprocess.getoutput(cmd_sub)	

	cmd_int = bdtArg + " intersect -header -wa -wb"
	if sortArg:
		cmd_int += " -sorted"
	cmd_int += " -a " + outputArg + "__tmp__.sorted.vcf"
	cmd_int += " -b " + rmArg
	cmd_int += " > " + outputArg + "__tmp_intersected__.vcf"
	print (cmd_int, "\n")
	subprocess.getoutput(cmd_int)

	writeLines = list()
	with open(outputArg + "__tmp_intersected__.vcf", "r") as vcfLines:
		for vcfLine in vcfLines:
			if vcfLine.startswith("#"):
				continue
			vcfLine = vcfLine.rstrip().split("\t")

			ref = vcfLine[3]
			alt = vcfLine[4]
			refdb = vcfLine[13]
			altdb = vcfLine[14]

			if len(list(refdb)) == 1:
				if ref == refdb:
					if len(list(altdb)) == 1:
						if alt != altdb:
							writeLine = "\t".join(vcfLine[:10]) + "\n"
							writeLines.append(writeLine)

					elif len(altdb.split(",")) > 1:
						altdb = altdb.split(",")
						if alt not in altdb:
							writeLine = "\t".join(vcfLine[:10]) + "\n"
							writeLines.append(writeLine)
	
	outwrite = open(outputArg + "__tmp_intersected_save__.vcf", "w")
	outwrite.write("".join(writeLines))
	outwrite.close()

	cmd_merge = "cat " + outputArg + "__tmp_subtracted__.vcf" + " " + outputArg + "__tmp_intersected_save__.vcf" + " > " + outputArg + "__tmp_merged__.vcf"
	print (cmd_merge, "\n")
	subprocess.getoutput(cmd_merge)

	resultFile = vcfArg.split("/")[-1].rstrip("vcf").rstrip("flt")
	cmd_sort = 'grep "^#" ' + outputArg + "__tmp_merged__.vcf" + " > " + outputArg + resultFile + "flt.vcf"
	cmd_sort += ' && grep -v "^#" ' + outputArg + "__tmp_merged__.vcf" + " | sort -V -k1,1 -k2,2n >> "  + outputArg + resultFile + "flt.vcf"
	print (cmd_sort, "\n")
	subprocess.getoutput(cmd_sort)

	cmd_rm = remove_temp(outputArg, "__tmp")
	print (cmd_rm)
	subprocess.getoutput(cmd_rm)

	cmd_rm = remove_temp("", vcfArg)
	print (cmd_rm, "\n")
	subprocess.getoutput(cmd_rm)

	return outputArg + resultFile + "flt.vcf"


def flt_vaf(vcfArg, vafArg):
	openVcf = open(vcfArg, "r")
	vcfLines = openVcf.readlines()
	openVcf.close()
	
	writeLines = [x for x in vcfLines if x.startswith("#")]
	for vcfLine in vcfLines:
		if vcfLine.startswith("#"):	continue
		tmpLine = vcfLine.rstrip().split()
		info = tmpLine[-1]
		depth = int(info.split(":")[0])
		variant = int(info.split(":")[1].split(",")[1])
		vaf = variant/depth*100.0
		if vaf <= vafArg:
			writeLines.append(vcfLine)

	outwrite = open(vcfArg, "w")
	outwrite.write("".join(writeLines))
	outwrite.close()
	
	return vcfArg


def tmp_vcf(vcfArg):
	tmpVcfName = vcfArg.replace("flt.vcf", "vcf")

	cmd_tmpvcf = "mv " + vcfArg + " " + tmpVcfName
	print (cmd_tmpvcf, "\n")

	subprocess.getoutput(cmd_tmpvcf)
	return tmpVcfName 



