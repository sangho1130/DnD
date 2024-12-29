
# 2024-12-28; version 1.0

import os
import sys
import io
import argparse
import subprocess
import pysam
import pandas as pd
import numpy as np


# 1
def Convert(t, d):
	for a, b in t:
		d.setdefault(b, []).append(a)
	return d


def get_read_info(chunk, bam_path):
	samfile = pysam.AlignmentFile(bam_path, "rb")
	read_level = []
	bed_df =[]

	umiFlag = False
	for variant in chunk.index:
		REAL_POSITION= chunk.loc[variant,"position"]
		PYSAM_POSITION = int(REAL_POSITION) - 1
        
		variantAllele = chunk.loc[variant,"alt"]
		chrom = chunk.loc[variant,"chr"]

		for read in samfile.fetch(chrom, PYSAM_POSITION, REAL_POSITION):
			if read.is_duplicate == True:	continue
            		
			di = {}
			di = Convert(read.get_aligned_pairs(),di)

			if type(di[PYSAM_POSITION][0]) == int:
				read_seq = read.query_sequence
				BASE = read_seq[di[PYSAM_POSITION][0]]

				if BASE == variantAllele:
					name = read.query_name
					add_seq = "="*di[PYSAM_POSITION][0] + BASE + "="*len(read_seq[di[PYSAM_POSITION][0]+1: ])
					read_flag = read.flag
					read_bc = read.get_tag("CB")
					read_level.append([chunk.loc[variant,'change'], name, str(read_flag), add_seq, read_bc])
	
	if umiFlag:
		readCols =["variantName", "readName", "readFlag", "readSeq", "Barcode", "UMI"]
	else:
		readCols =["variantName", "readName", "readFlag", "readSeq", "Barcode"]
	read_level = pd.DataFrame(read_level, columns = readCols)
	read_level = read_level.values.tolist()
	
	return (read_level)
	

def select_redun_reads(read_level):
	mateDict = dict()
	editDict = dict()	
	for read_info in read_level:
		if read_info[1] not in mateDict:
			mateDict[read_info[1]] = dict()
			mateDict[read_info[1]][int(read_info[2])] = 1
		elif int(read_info[2]) not in mateDict[read_info[1]]:
			mateDict[read_info[1]][int(read_info[2])] = 1
		else:
			mateDict[read_info[1]][int(read_info[2])] += 1
	
		if read_info[0] not in editDict:
			editDict[read_info[0]] = dict()
			editDict[read_info[0]][read_info[1]] = 1
		elif read_info[1] not in editDict[read_info[0]]:
			editDict[read_info[0]][read_info[1]] = 1
		else:
			editDict[read_info[0]][read_info[1]] += 1

	editLines = list()
	for read_info in read_level:
		localDict = editDict[read_info[0]]
		redun_reads = [x[0] for x in localDict.items() if x[1] != 1]	

		if read_info[1] in redun_reads:
			mate_info = mateDict[read_info[1]]
			max_count = max(mate_info.values())
			if max_count == 1:
				max_flag = min(mate_info.keys())
			else:
				max_flag = sorted([x[0] for x in mate_info.items() if x[1] == max(mate_info.values())])[0]

			if int(read_info[2]) == max_flag:
				editLines.append(read_info)
		else:
			editLines.append(read_info)

	return (editLines)


def add_editColumn(readLines):
	seqDict = dict()
	for readLine in readLines:

		if readLine[1] not in seqDict:
			seqDict[readLine[1]] = dict()
			seqDict[readLine[1]][ int(readLine[2]) ] = readLine[3]
		elif int(readLine[2]) not in seqDict[readLine[1]]:
			seqDict[readLine[1]][ int(readLine[2]) ] = readLine[3]
		else:
			stored = list(seqDict[readLine[1]][ int(readLine[2]) ])
			newedit = list(readLine[3])
			updated = ""
			for i in range(len(stored)):
				if stored[i] == newedit[i]:	updated += stored[i]
				elif stored[i] != "=":	updated += stored[i]
				elif newedit[i] != "=":	updated += newedit[i]
			seqDict[readLine[1]][ int(readLine[2]) ] = updated

	editLines = list()
	for readLine in readLines:
		seqs = seqDict[readLine[1]]
		if len(seqs.keys()) == 1:
			readLine.append(seqDict[readLine[1]][ int(readLine[2]) ])
		else:
			addseq = list()
			for flag in sorted(seqs.keys()):
				addseq.append( seqs[flag] )
			addseq = "|".join(addseq)
			readLine.append(addseq)
		editLines.append(readLine)

	return (editLines)


def flt_bam_qnames(smtArg, bamArg, qnameArg):
	suffix = qnameArg.split("/")[-1].split(".")[0]
	cmd_flt = smtArg + " view -b"
	cmd_flt += " -N " + qnameArg
	cmd_flt += " -o " + bamArg[:-3] + suffix + "_edits.bam"
	cmd_flt += " " + bamArg
	subprocess.getoutput(cmd_flt)

	cmd_idx = smtArg + " index " + bamArg[:-3] + suffix + "_edits.bam"
	subprocess.getoutput(cmd_idx)


def make_edit_info(smtArg, bgzipArg, tabixArg, vcfArg, bamArg, outputArg, refArg, altArg, fltbamArg):
	edits = pd.read_table(vcfArg, header = None, comment='#', usecols = [0,1,3,4])
	edits.columns = ["chr", "position", "ref", "alt"]
	edits["change"] = edits["chr"] + ">" + edits["position"].astype(str) + ">" + edits["ref"] + ">" + edits["alt"]
	edits.index = edits.loc[:, "change"]

	revcompDict = {"A":"T", "T":"A", "G":"C", "C":"G"}
	edits = edits[ (edits["ref"] == refArg) & (edits["alt"] == altArg) | (edits["ref"] == revcompDict[refArg]) & (edits["alt"] == revcompDict[altArg])]

	readNameLines = get_read_info(edits, bamArg)
	writeLines = select_redun_reads(readNameLines)
	writeLines = add_editColumn(writeLines)

	fragDict = dict()
	for writeLine in writeLines:
		checkVar = writeLine[0].split(">")[0] + "\t" + str(int(writeLine[0].split(">")[1])-1) + "\t" + writeLine[0].split(">")[1] + "\t" + writeLine[4]
		if checkVar not in fragDict:
			fragDict[checkVar] = 1
		else:
			fragDict[checkVar] += 1

	outputArg += "dndedits.txt"

	print ("writing fragments.tsv")	
	outwrite = open(outputArg, "w")
	outwrite_f = open("/".join(outputArg.split("/")[:-1]) + "/fragments.tsv", "w")
	prevLineCheck = ""
	for writeLine in writeLines:
		outwrite.write("\t".join(writeLine) + "\n")
		checkVar = writeLine[0].split(">")[0] + "\t" + str(int(writeLine[0].split(">")[1])-1) + "\t" + writeLine[0].split(">")[1] + "\t" + writeLine[4]
		fragmentLine = checkVar + "\t" + str(fragDict[checkVar])
		if fragmentLine == prevLineCheck:
			continue
		prevLineCheck = fragmentLine
		outwrite_f.write(fragmentLine + "\n")

	outwrite.close()
	outwrite_f.close()

	cmd_gz = bgzipArg + " " + "/".join(outputArg.split("/")[:-1]) + "/fragments.tsv"
	print (cmd_gz)
	subprocess.getoutput(cmd_gz)
	cmd_idx = tabixArg + " --preset=bed " + "/".join(outputArg.split("/")[:-1]) + "/fragments.tsv.gz"
	print (cmd_idx)
	subprocess.getoutput(cmd_idx)

	qnames = sorted(list(set( [x[1] for x in writeLines] )))
	outwrite = open(outputArg[:-3] + "qnames.txt", "w")
	outwrite.write("\n".join(qnames) + "\n")
	outwrite.close()
	
	if fltbamArg:
		print ("filtering bam")
		flt_bam_qnames(smtArg, bamArg, outputArg[:-3] + "qnames.txt")

	return (outputArg)


# 2
def make_fragments(bdtArg, peakArg, editArg, outputArg):
	print ("summarizing edits")
	openPeak = open(peakArg, "r")
	peakLines = openPeak.readlines()
	openPeak.close()

	peakChrDict = dict()
	for peakLine in peakLines:
		peak_chr = peakLine.split()[0]
		if peak_chr not in  peakChrDict:
			peakChrDict[peak_chr] = list()
		peakChrDict[peak_chr].append(peakLine)

	openEdit = open(editArg, "r")
	editLines = openEdit.readlines()
	openEdit.close()

	editChrDict = dict()
	for editLine in editLines:
		edit_chr = editLine.split()[0].split(">")[0]
		if edit_chr not in  editChrDict:
			editChrDict[edit_chr] = list()
		editChrDict[edit_chr].append(editLine)

	chrmOrder = ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10',
                        'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
                        'chr21', 'chr22', 'chrX', 'chrY']
	chrmOrder = [x for x in chrmOrder if x in peakChrDict.keys()]
	chrmOrder = [x for x in chrmOrder if x in editChrDict.keys()]

	countDict = dict()
	bcCountDict = dict() ###
	for chrm in chrmOrder:
		peakLines = peakChrDict[chrm] 
		editLines = editChrDict[chrm]

		for peakLine in peakLines:
			peakLine = peakLine.rstrip().split()

			peakLine_chr = peakLine[0]
			peakLine_start = int(peakLine[1])
			peakLine_end = int(peakLine[2])

			i = 0
			for editLine in editLines:
				editLine = editLine.rstrip().split()

				edit_chr = editLine[0].split(">")[0]
				edit_pos = int(editLine[0].split(">")[1])
				edit_qname = editLine[1]
				edit_bc = editLine[4]

				if peakLine_chr == edit_chr:
					if peakLine_start <= edit_pos <= peakLine_end:
						if "_".join(peakLine[:3]) not in countDict:
							countDict["_".join(peakLine[:3])] = dict()
							countDict["_".join(peakLine[:3])][edit_qname + "_" + edit_bc] = 1
						elif edit_qname + "_" + edit_bc not in countDict["_".join(peakLine[:3])]:
							countDict["_".join(peakLine[:3])][edit_qname + "_" + edit_bc] = 1
						else:
							countDict["_".join(peakLine[:3])][edit_qname + "_" + edit_bc] += 1
						
						if edit_bc not in bcCountDict:
							bcCountDict[edit_bc] = 1
						else:	bcCountDict[edit_bc] += 1
						i += 1
					else:
						break
				else:
					break
	
			editLines = editLines[i:]

	tmpDict = dict()
	for coordKey in sorted(countDict.keys()):
		qname_bcs = sorted(countDict[coordKey].keys())
		for qname_bc in qname_bcs:
			tmpKey = coordKey.replace("_", "\t") + "\t" + qname_bc.split("_")[1]
			if tmpKey not in tmpDict:
				tmpDict[tmpKey] = countDict[coordKey][qname_bc]
			else:
				tmpDict[tmpKey] += countDict[coordKey][qname_bc]
	
	writeLines = list()
	for coordKey in sorted(countDict.keys()):
		qname_bcs = sorted(countDict[coordKey].keys())
		for qname_bc in qname_bcs:
			tmpKey = coordKey.replace("_", "\t") + "\t" + qname_bc.split("_")[1]
			if tmpKey in tmpDict:
				writeLine = tmpKey + "\t" + str(tmpDict[tmpKey]) + "\n"
				writeLines.append(writeLine)
				tmpDict.pop(tmpKey)	

	outputArg += "dndedits.fragments.txt"
	outwrite = open(outputArg + ".unsorted", "w")
	outwrite.write("".join(writeLines))
	outwrite.close()

	cmd_bdt = bdtArg + " sort -i " + outputArg + ".unsorted"
	cmd_bdt += " > " + outputArg
	subprocess.getoutput(cmd_bdt)

	cmd_rm = "rm " + outputArg + ".unsorted"
	subprocess.getoutput(cmd_rm)

	outwrite = open(outputArg[:-3] + "editCount.txt", "w")
	for bc in sorted(bcCountDict.keys()):
		writeLine = bc + "\t" + str(bcCountDict[bc]) + "\n"
		outwrite.write(writeLine)
	outwrite.close()

	return (outputArg)


# 3
def make_mtx(fragArg, outputArg):
	print ("writing mtx")
	outputArg += "mtx/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	openFrag = open(fragArg, "r")
	fragLines = openFrag.readlines()
	openFrag.close()

	peakLines = list()
	prev_peakLine = ""

	barcodes = set()

	countDict = dict()
	allcounts = 0
	for fragLine in fragLines:
		fragLine = fragLine.rstrip().split()
		
		barcodes.add(fragLine[3])		

		peakLine = "\t".join(fragLine[:3])
		if peakLine != prev_peakLine:
			peakLines.append(peakLine)
			prev_peakLine = peakLine

		if fragLine[3] not in countDict:
			countDict[fragLine[3]] = dict()
		countDict[fragLine[3]][peakLine] = fragLine[-1]
		allcounts += int(fragLine[-1])

	barcodes = sorted( list(barcodes) )

	outwrite = open(outputArg + "peaks.bed", "w")
	outwrite.write("\n".join(peakLines) + "\n")
	outwrite.close()

	outwrite = open(outputArg + "barcodes.tsv", "w")
	outwrite.write("\n".join(barcodes) + "\n")
	outwrite.close()

	writeLines = list()
	writeLines.append("%%MatrixMarket matrix coordinate integer general\n")
	writeLines.append("%dnd_custom\n")
	writeLines.append(str(len(peakLines)) + "\t" + str(len(barcodes)) + "\t" + str(allcounts) + "\n")

	for b in range(len(barcodes)):
		bc = barcodes[b]
		bcCounts = countDict[bc]
		
		for p in range(len(peakLines)):
			peak = peakLines[p]
			if peak in bcCounts:
				writeLine = str(p+1) + "\t" + str(b+1) + "\t" + bcCounts[peak] + "\n"
				writeLines.append(writeLine)
	outwrite = open(outputArg + "matrix.mtx", "w")
	outwrite.write("".join(writeLines))
	outwrite.close()



def main(args):
	smt = subprocess.getoutput("which samtools")
	bgzip = subprocess.getoutput("which bgzip")
	tabix = subprocess.getoutput("which tabix")
	bdt = subprocess.getoutput("which bedtools")

	if args.Output[-1] != "/":
		args.Output += "/"

	# 1
	editInfo = make_edit_info(smt, bgzip, tabix, args.Vcf, args.Bam, args.Output, args.ref, args.alt, args.fltbam)
	# 2
	fragInfo = make_fragments(bdt, args.Peak, editInfo, args.Output)
	# 3
	make_mtx(fragInfo, args.Output)


if __name__ == '__main__':
	parser = argparse.ArgumentParser('')
	parser.add_argument('-v', '--Vcf', help = 'input vcf file', required = True)
	parser.add_argument('--ref', help = "(*optional) reference allele, rev.comp. is automatically considered; default is C", required = False, default = "C")
	parser.add_argument('--alt', help = "(*optional) altered allele, rev.comp. is automatically considered; default is T", required = False, default = "T")

	parser.add_argument('-p', '--Peak', help = 'input peak file', required = True)
	parser.add_argument('-b', '--Bam', help = 'input bam file', required = True)
	parser.add_argument('--flt-bam', dest = "fltbam", help = "(*optional) filter <-b/--Bam> using variants; default is NOT FILTERING", required = False, default = False, action = 'store_true')

	parser.add_argument('-o', '--Output', help = 'output directory path', required = True)

	args = parser.parse_args()
	main(args)



