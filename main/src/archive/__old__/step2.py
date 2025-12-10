
import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process


### re-write to multi-core compatible version
def mpileup_make(bcfArg, threadArg, outputArg, bamArgs, fastaArg):
	outputArg += "step2_mpileup/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	cmd = bcfArg + " mpileup"
	cmd += " -I"
	cmd += " -a FORMAT/AD,FORMAT/DP,INFO/AD"
	cmd += " --no-BAQ --min-MQ 1"
	cmd += " --max-depth 8000"
	cmd += " --threads " + threadArg
	cmd += " -Ov -f " + fastaArg 

	mpileups = list()
	for bamArg in bamArgs:
		sampleName = bamArg.split("/")[-2]
		if not os.path.exists(outputArg + sampleName):
			os.mkdir(outputArg + sampleName)

		cmd_pu = cmd + " " + bamArg
		cmd_pu += " -o " + outputArg + sampleName + "/" + bamArg.split("/")[-1].replace("bam", "pileup")
		mpileups.append(outputArg + sampleName + "/" + bamArg.split("/")[-1].replace("bam", "pileup"))
		print (cmd_pu, "\n")
		subprocess.getoutput(cmd_pu)

	return mpileups


def mpileup_vcf(pileupArg, countArg, altArg):
	writeLines = list()

	with open(pileupArg) as vcfLines:
		for vcfLine in vcfLines:
			if vcfLine.startswith("#"):
				writeLines.append(vcfLine)
				continue

			vcfLine_s = vcfLine.rstrip().split("\t")
			if vcfLine_s[0] == "chrM":	continue
			if not vcfLine_s[0].startswith("chr"):	continue
			if vcfLine_s[3] == "N":	continue

			alt = vcfLine_s[4].split(",")
			if len(alt) == 1:	continue
			if alt[0] == "N":	continue

			summ = vcfLine_s[-1].split(":")
			summ_dp = int(summ[1])
			if summ_dp < countArg:  continue

			summ_counts = [int(x) for x in summ[2].split(",")]
			info = vcfLine_s[7].split(";")
			info_ad = info[1].lstrip("AD=").split(",")
			info_qs = info[3].lstrip("QS=").split(",")
			info_qs = [float(x) for x in info_qs]

			alt = [a for a in alt if a != "<*>"]

			for i in range(len(alt)):
				if summ_counts[i+1] < altArg:	continue
				
				writeLine = vcfLine_s[:4]
				writeLine.append(alt[i])

				writeLine += vcfLine_s[5:7]
				
				a_info = list()
				a_info.append( info[0] )
				a_info.append( "AD=" + ",".join( [info_ad[0], info_ad[i+1]] ) )
				a_info.append( "QS=" + ",".join( [str(info_qs[0]), str(info_qs[i+1])] ) )
				a_info = ";".join(a_info)

				writeLine.append(a_info)
				
				writeLine.append( "DP:AD" )
				writeLine.append( summ[1] + ":" + str(summ_counts[0]) + "," + str(summ_counts[i+1]) )
				
				writeLine = "\t".join(writeLine) + "\n"
				writeLines.append(writeLine)

	outputFile = pileupArg.replace(".pileup", ".vcf")
	outwrite = open(outputFile, "w")
	outwrite.write("".join(writeLines))
	outwrite.close()

	return outputFile



