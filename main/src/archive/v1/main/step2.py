
import argparse
import os
import sys
import subprocess
import gzip
import random
import pysam
from Bio import SeqIO
from multiprocessing import Pool


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


def pysam_pileup(bamArg, chrArg, countArg, altCountArg, l_trimArg, r_trimArg, gDict):
	samfile = pysam.AlignmentFile(bamArg, "rb")

	varLines = list()
	for pileupcolumn in samfile.pileup(chrArg):
		pos = pileupcolumn.pos +1
		dp = pileupcolumn.n
		if dp < countArg: continue

		ref_allele = gDict[chrArg].seq[pileupcolumn.pos]
		if ref_allele == "N": continue

		ntDict = {"A":0, "T":0, "G":0, "C":0}
		for pileupread in pileupcolumn.pileups:
			pos_in_read = pileupread.query_position
			if pos_in_read == None: continue # deletion

			nt_in_read = pileupread.alignment.query_sequence[pos_in_read]
			if nt_in_read == "N": continue
		
			len_of_read = len(pileupread.alignment.query_sequence)

			if nt_in_read != ref_allele:
				if pos_in_read > l_trimArg and pos_in_read < len_of_read - r_trimArg:
					ntDict[nt_in_read] += 1
				else:	ntDict[ref_allele] += 1
			else:	ntDict[ref_allele] += 1

		if sum(ntDict.values()) < countArg:	continue

		ref_allele_count = str(ntDict[ref_allele])
		alt_alleles = [alt for alt in ntDict.keys() if alt != ref_allele and ntDict[alt] >= altCountArg]
		for alt_allele in alt_alleles:
			if ntDict[alt_allele] == 0:	continue
			pos_info = "DP=" + str(dp) + ";" + "AD=" + str(ntDict[ref_allele]) + "," + str(ntDict[alt_allele])
			alt_info = str(ntDict[ref_allele] + ntDict[alt_allele]) + ":" + str(ntDict[ref_allele]) + "," + str(ntDict[alt_allele])

			writeLine = [chrArg, str(pos), ".", ref_allele, alt_allele, "0", ".", pos_info, "DP:AD", alt_info]
			writeLine = "\t".join(writeLine) + "\n"
			varLines.append(writeLine)

	samfile.close()
	return (varLines)


def pysam_pileup_mp(args):
	bamArg, chrm, countArg, altArg, leftTrimArg, rightTrimArg, fastaArg = args
	record_dict = SeqIO.index(fastaArg, "fasta")
	puLines_chrm = pysam_pileup(bamArg, chrm, countArg, altArg, leftTrimArg, rightTrimArg, record_dict)
	return (puLines_chrm)


def pu_vcf_make(bamArgs, threadArg, outputArg, fastaArg, chromArg, countArg, altArg, leftTrimArg, rightTrimArg):
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)
	outputArg += "step3_fltvcf/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	record_dict = SeqIO.index(fastaArg, "fasta")
	chrms = list(record_dict.keys())
	if not chromArg:
		chrms = [x for x in chrms if x.startswith("chr") and x != "chrM"]
	chrms.sort()

	resvcfs = list()
	for bamArg in bamArgs:
		sampleName = bamArg.split("/")[-2]
		if not os.path.exists(outputArg + sampleName):
			os.mkdir(outputArg + sampleName)
		resvcf = outputArg + sampleName + "/" + bamArg.split("/")[-1].replace("bam", "vcf")
		
		writeLines = list()
		writeLines.append("##fileformat=VCFv4.2\n")
		writeLines.append("##reference=file://" + fastaArg + "\n")

		for chrm in chrms:
			writeLine = "##contig=<ID=" + chrm + ",length=" + str(len(record_dict[chrm])) + ">\n"
			writeLines.append(writeLine)
		writeLines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tPILEUP\n")

		args_list = [ (bamArg, chrm, countArg, altArg, leftTrimArg, rightTrimArg, fastaArg) for chrm in chrms ]
		p = Pool(threadArg)
		try:
			print ("running")
			puLines = p.map(pysam_pileup_mp, args_list)
		finally:
			print ("pooling")
			p.close()
			p.join()
		print ("sum")
		puLines = sum(puLines, [])
		print ("append")
		writeLines += puLines

		outwrite = open(resvcf, "w")
		outwrite.write("".join(writeLines))
		outwrite.close()
		resvcfs.append(resvcf)

	return (resvcfs)



