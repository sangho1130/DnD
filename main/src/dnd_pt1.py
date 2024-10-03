
import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process
import time


###################################
### Step 1. preprocess bam file ###
##################################

def picard_markDup(pcdArg, dirArg, outputArg, bamArg):
	cmd = pcdArg + " MarkDuplicates"
	cmd += ' I=' + dirArg + bamArg
	cmd += ' O=' + outputArg + "__tmp_rmDup__" + bamArg.lower()
	cmd += ' M=' + outputArg + "step1_picard_deduplication.txt" 
	cmd += ' REMOVE_DUPLICATES=true'
	return cmd, "__tmp_rmDup__" + bamArg.lower()

def samtools_index(smtArg, outputArg, bamArg):
	if outputArg == "":
		cmd = smtArg + ' index ' + bamArg
	else:
		cmd = smtArg + ' index ' + outputArg + bamArg
	return cmd

def samtools_flt(smtArg, outputArg, bamArg, threadArg, mapqArg, chromArg, otherArg, seArg):
	cmd = smtArg + " view -b"
	cmd += ' --threads ' + threadArg
	cmd += ' -q ' + mapqArg # 20 by default
	if not seArg:
		cmd += ' -f 0x2' # read mapped in proper pair
	cmd += ' -F 256' # leave only primary alignment
	if not chromArg:
		cmd += ' -L ' + 'dnd_acc/GRCh38.p14.genome.chrs.bed' # Leave only chromosomes
	if otherArg:
		cmd += ' ' + otherArg
	cmd += ' ' + outputArg + bamArg
	cmd += ' > ' + outputArg + bamArg.replace('__tmp_rmDup__', '')
	return cmd, bamArg.replace('__tmp_rmDup__', '')

def remove_temp(pathArg, patternArg):
	if pathArg == "":
		pathArg = patternArg
	elif pathArg[-1] != "/":
		pathArg += "/" + patternArg
	else:
		pathArg += patternArg
	cmd = "rm " + pathArg + "*"
	return cmd

def flt_bam(smtArg, pcdArg, dirArg, outputArg, threadArg, mapqArg, chromArg, otherArg, seArg, startArg, endArg):
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)
	outputArg += "step1_preprocess/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	bamDirs = os.listdir(dirArg)
	bamDirs.sort()
	
	if startArg:
		if endArg:
			bamDirs = bamDirs[startArg:endArg]
		else:	bamDirs = bamDirs[startArg:]
	else:
		if endArg:
			bamDirs = bamDirs[:endArg]
		else:	pass

	resbams = list()
	for bamDir in bamDirs:
		bamFile = [x for x in os.listdir(dirArg + bamDir) if x.endswith(".bam")][0]
		if not os.path.exists(outputArg + bamDir):
			os.mkdir(outputArg + bamDir)
		cmd_markDup, resBam = picard_markDup(pcdArg, dirArg + bamDir + "/", outputArg + bamDir + "/", bamFile)
		print (cmd_markDup, "\n")
		subprocess.getoutput(cmd_markDup)
		
		cmd_smtFlt, resBam = samtools_flt(smtArg, outputArg + bamDir + "/", resBam, threadArg, mapqArg, chromArg, otherArg, seArg)
		print (cmd_smtFlt, "\n")
		subprocess.getoutput(cmd_smtFlt)

		cmd_smtIdx = samtools_index(smtArg, outputArg + bamDir + "/", resBam)
		print (cmd_smtIdx, "\n")
		subprocess.getoutput(cmd_smtIdx)

		rm_cmd = remove_temp(outputArg + bamDir, "__tmp_rmDup__")
		print (rm_cmd, "\n")
		subprocess.getoutput(rm_cmd)

		resbams.append(outputArg + bamDir + "/" + resBam)
	return resbams


######################
### Step 2. pileup ###
######################

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


#############################
### Step 3. filtering vcf ###
#############################

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


#############################
### Step 4. extract edits ###
#############################

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


############################
### Step 5. Peak calling ###
############################

def macs2_callpeak(mcsArg, bdtArg, bamArg, outputArg, gsizeArg, optArg, bklist, passBkArg):
	outputArg += "step5_peaks/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)
	
	sampleName = bamArg.split("/")[-2]
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

	cmdmc2 = cmd + " -t " + bamArg
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
		mc2dir = macs2_callpeak(mc2, bdt, bam, args.Output, args.gsize, args.opt, args.blacklist, args.pass_bklist)
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
	parser.add_argument('--fasta', help = '[Step 2] genome fasta (indexed) used in alignment; e.g. cellranger/fasta/genome.fa', required = True)

	parser.add_argument('--gnomad', help = '[Step 3] (*optional) path to gnomAD vcf file', default = "dnd_acc/somatic-hg38-af-only-gnomad.hg38.chrs.vcf.gz", required = False)
	parser.add_argument('--pass-gnomad', dest = 'pass_germline', help = '[Step 3] (*optional) do not run gnomAD filering', required = False, default = False, action = 'store_true')
	parser.add_argument('--custom', help = '[Step 3] (*optional) custom vcf file(s) to filter, multiple files are accepted', nargs = "+", required = False)
	parser.add_argument('--vaf', help = '[Step 3] (*optional) filter mutations frequent than this value; e.g. 10 for 10perc; default is 10', type = int, default = 10, required = False)

	parser.add_argument('--snv', help = '[Step 4] (*optoinal) SNV patterns to search; e.g. --snv "C>T,G>A, G>C"; default is "C>T,G>A"', default = "C>T,G>A", required = False)

	parser.add_argument('--gsize', help = '[Step 5] (*optional) effective genome size for macs2 callpeak; default is hs (homo sapiens)', default = 'hs', required = False)
	parser.add_argument('--opt', help = '[Step 5] (*optional) other parameters for macs2 callpeak', type = str, required = False)
	parser.add_argument('--blacklist', help = '[Step 5] (*optional) blacklist file; default is hg38-blacklist.v2.bed', required = False, default = "dnd_acc/hg38-blacklist.v2.bed")
	parser.add_argument('--pass-bklist', dest = 'pass_bklist', help = '[Step 5] (*optional) do not run blacklist filtering', required = False, default = False, action = 'store_true')

	args = parser.parse_args()
	main(args)


