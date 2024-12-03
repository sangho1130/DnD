
import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process
from main.dnd_acc import *


def flt_snvs(vcfArg, outputArg, snvArg):
	outputArg += "step4_snvs/"
	if not os.path.exists(outputArg):
                os.mkdir(outputArg)

	sampleName = vcfArg.split("/")[-2]
	outputArg += sampleName + "/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	outVcfArg = outputArg + vcfArg.split("/")[-1]
	outVcfArg = outVcfArg.rstrip("\\.vcf.gz")
	if not snvArg:
		outVcfArg += ".SNVs.vcf"
	else:
		outVcfArg += "." + snvArg.replace(" ", "").replace(",", "_").replace(">", "") + ".SNVs.vcf"

	snvDict = dict()
	if not snvArg:
                snvDict["A"] = ["T", "G", "C"]
                snvDict["T"] = ["A", "G", "C"]
                snvDict["G"] = ["A", "T", "C"]
                snvDict["C"] = ["A", "T", "G"]
	else:
                snvArg = snvArg.strip('"')
                snvArg = snvArg.split(",")
                snvArg = [x.strip(" ") for x in snvArg]
                for snv in snvArg:
                        ref = snv.split(">")[0]
                        alt = snv.split(">")[1]
                        if not ref in snvDict:
                                snvDict[ref] = [alt]
                        else:
                                snvDict[ref].append(alt)
	print ("Searching for {'Ref': ['Alt'], ... }:", snvDict, "\n")
	
	if vcfArg.endswith(".vcf.gz"):
                openVcf = gzip.open(vcfArg, "rb")
	else:   openVcf = open(vcfArg, "r")
	vcfLines = openVcf.readlines()
	openVcf.close()

	writeLines = list()
	for vcfLine in vcfLines:
                if vcfLine.startswith("#"):
                        writeLines.append(vcfLine)
                        continue

                refAll = vcfLine.rstrip().split()[3]
                altAll = vcfLine.rstrip().split()[4]

                if refAll in snvDict.keys():
                        if altAll in snvDict[refAll]:
                                writeLines.append(vcfLine)

	outwrite = open(outVcfArg, "w")
	outwrite.write(''.join(writeLines))
	outwrite.close()

	return outVcfArg


def flt_alignment(bdpArg, smtArg, bamArg, vcfArg, outputArg):
	outputArg += "step4_snvs/"
	sampleName = bamArg.split("/")[-2]
	outputArg += sampleName + "/"

	outbedArg = outputArg + bamArg.split("/")[-1].replace(".bam", ".snvs.bed")
	outbamArg = outputArg + bamArg.split("/")[-1].replace(".bam", ".snvs.bam")

	BAM_FILE = 'BAM_FILE="' + bamArg + '"'
	VCF_FILE = 'VCF_FILE="' + vcfArg + '"'
	OUTPUT_FILE = 'OUTPUT_FILE="' + outbedArg + '"'

	cmd_bedops = bdpArg
	cmd_bedops += " -e 1"
	cmd_bedops += ' <(bam2bed < "$BAM_FILE")'
	cmd_bedops += ' <(vcf2bed < "$VCF_FILE")'
	cmd_bedops += ' > "$OUTPUT_FILE"'
	print (BAM_FILE)
	print (VCF_FILE)
	print (OUTPUT_FILE)
	print (cmd_bedops, "\n")

	randsig = random.random()
	if len(bamArg.split("/")) != 1:
		randFname = outputArg + "__temp__" + str(randsig) + ".sh"
	else:
		pwd = os.getcwd()
		randFname = pwd + "/__temp__" + str(randsig) + ".sh"

	openFile = open(randFname, "w")
	openFile.write("#!/bin/bash -l\n\n")
	openFile.write("mamba activate scatacseq\n")
	openFile.write(BAM_FILE + "\n")
	openFile.write(VCF_FILE + "\n")
	openFile.write(OUTPUT_FILE + "\n")
	openFile.write(cmd_bedops + "\n")
	openFile.close()

	subprocess.getoutput("chmod +x " + randFname)
	subprocess.getoutput(randFname)
	subprocess.getoutput("rm " + randFname)

	cmd_smt = smtArg
	cmd_smtvw = cmd_smt + " view -b -h"
	cmd_smtvw += " -N " + outbedArg 
	cmd_smtvw += " " + bamArg
	cmd_smtvw += " > " + outbamArg
	print (cmd_smtvw, "\n")
	subprocess.getoutput(cmd_smtvw)

	cmd_smtid = samtools_index(smtArg, "", outbamArg)
	print (cmd_smtid, "\n")
	subprocess.getoutput(cmd_smtid)

	cmd_rm = remove_temp("", outbedArg)
	subprocess.getoutput(cmd_rm)

	return outbamArg



