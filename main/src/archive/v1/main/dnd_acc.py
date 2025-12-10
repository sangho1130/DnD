
import argparse
import os
import sys
import subprocess
import gzip
import random


def remove_temp(pathArg, patternArg):
	if pathArg == "":
		pathArg = patternArg
	elif pathArg[-1] != "/":
		pathArg += "/" + patternArg
	else:
		pathArg += patternArg
	cmd = "rm " + pathArg + "*"
	return cmd


def samtools_index(smtArg, outputArg, bamArg):
	if outputArg == "":
		cmd = smtArg + ' index ' + bamArg
	else:
		cmd = smtArg + ' index ' + outputArg + bamArg
	return cmd


def fltRedunPeaks(peakArg):
	openPeak = open(peakArg, "r")
	peakLines = openPeak.readlines()
	openPeak.close()

	uniquePeaks = list(set(["_".join(x.split()[:3]) for x in peakLines]))
	if len(uniquePeaks) == len(peakLines): return
	n = 1
	used_peaks = list()
	writeLines = list()
	for peakLine in peakLines:
		peakLine = peakLine.rstrip().split()
		tmpInfo = "_".join(peakLine[:3])
		if tmpInfo not in used_peaks:
			used_peaks.append(tmpInfo)
			peakLine[3] = "peak_" + str(n)
			writeLines.append("\t".join(peakLine) + "\n")
			n += 1
	outwrite = open(peakArg, "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


def fltRedunVcfs(vcfArg, outputArg):
	openVcf = open(vcfArg, "r")
	vcfLines = openVcf.readlines()
	openVcf.close()

	if not outputArg:	outputArg = vcfArg[:-3] + "redunFlt.vcf"

	outwrite = open(outputArg, "w")
	redunCheck = list()
	for vcfLine in vcfLines:
		tmp = vcfLine.split()
		tmp_check = "_".join(tmp[0:2] + tmp[3:5])
		if tmp_check not in redunCheck:
			redunCheck.append(tmp_check)
			outwrite.write(vcfLine)
	outwrite.close()


def dnd_vcf(vcfArg, varArgs):
	outVcfArg = vcfArg.replace("SNVs.vcf", "_".join(varArgs) + ".SNVs.vcf")

	snvDict = dict()
	for varArg in varArgs:
		if varArg[0] not in snvDict:
			snvDict[ varArg[0] ] = [ varArg[1] ]
		else:	snvDict[ varArg[0] ].append(varArg[1])
	
	openVcf = open(vcfArg, "r")
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


def non_Dnd_vcf(vcfArg, varArgs):
	openFile = open(vcfArg, "r")
	fileLines = openFile.readlines()
	openFile.close()

	writeLines = list()
	for fileLine in fileLines:
		if fileLine.startswith("#"):
			writeLines.append(fileLine)
		elif fileLine.split()[3] + fileLine.split()[4] not in varArgs:
			writeLines.append(fileLine)

	outwrite = open(vcfArg.replace("vcf", "nonDnD.vcf"), "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


def merge_beds(motif_peak_files):
	openMp = open(motif_peak_files[0], 'r')
	mpLines = openMp.readlines()
	openMp.close()

	header = mpLines[0]
	mpLines = mpLines[1:]
	peaks = list(set([x.split()[0] for x in mpLines]))
	peaks.sort()

	for motif_peak_file in motif_peak_files[1:]:
		openFile = open(motif_peak_file, "r")
		fileLines = openFile.readlines()
		fileLines = fileLines[1:]
		openFile.close()

		for fileLine in fileLines:
			tmpLine = fileLine.rstrip().split('\t')
			if tmpLine[0] in peaks:	continue
			else:
				peaks.append(tmpLine[0])
				mpLines.append(fileLine)
	
	outwrite = open(motif_peak_files[0].replace("_1.txt", ".txt"), "w")
	outwrite.write(header)
	outwrite.write("".join(mpLines))
	outwrite.close()


def countOnlySnvs(vcfArg, peakArg, outputArg):
	if not outputArg:
		if vcfArg.endswith(".vcf.gz"):
			outputArg = vcfArg.rstrip("vcf.gz")
		else:	outputArg = vcfArg.rstrip("vcf")
	
	if vcfArg.endswith(".vcf.gz"):
		openVcf = gzip.open(vcfArg, "rb")
	else:	openVcf = open(vcfArg, "r")
	vcfLines = openVcf.readlines()
	vcfLines = [x for x in vcfLines if x[0] != "#"]
	openVcf.close()

	openPeak = open(peakArg, "r")
	peakLines = openPeak.readlines()
	openPeak.close()

	countDict = dict()
	ppeakDict = dict()
	revcompDict = {"T>A":"A>T", "T>C":"A>G", "T>G":"A>C",
			"G>A":"C>T", "G>T":"C>A", "G>C":"C>G",
			"A>T":"A>T", "A>G":"A>G", "A>C":"A>C",
			"C>T":"C>T", "C>A":"C>A", "C>G":"C>G"}

	for peakLine in peakLines:
		peakLine = peakLine.rstrip().split()
		chrm = peakLine[0]
		start = int(peakLine[1])
		end = int(peakLine[2])
		peakname = peakLine[3]

		peakrange = end-start
		peakKey = peakname + "_" + str(peakrange)
		ppeakDict[peakKey] = {"A>T":0, "A>G":0, "A>C":0, "C>A":0, "C>T":0, "C>G":0}

		for vcfLine in vcfLines:
			vcfLine = vcfLine.rstrip().split()

			chrm_vcf = vcfLine[0]
			coord_vcf = int(vcfLine[1])
			ref_vcf = vcfLine[3]
			alt_vcf = vcfLine[4]
			count_vcf = int(vcfLine[-1].split(",")[-1]) 

			if chrm != chrm_vcf:    continue
			elif start > coord_vcf: continue
			elif end < coord_vcf:   continue

			snpKey = ref_vcf + ">" + alt_vcf
			snpKey = revcompDict[snpKey]
			ppeakDict[peakKey][snpKey] += count_vcf

		for tmpKey in ppeakDict[peakKey].keys():
			if not tmpKey in countDict:	countDict[tmpKey] = 0.0
			if ppeakDict[peakKey][tmpKey] == 0:	continue
			countDict[tmpKey] += float(ppeakDict[peakKey][tmpKey])/len(peakLines)

	header = "Ref\tAlt\tSNVs_peak\n"
	outwrite = open(outputArg + ".countNormalized.txt", "w")
	outwrite.write(header)

	for snpKey in sorted(countDict.keys()):
		writeLine = "\t".join(snpKey.split(">")) + "\t" + str(countDict[snpKey]) + "\n"
		outwrite.write(writeLine)
	outwrite.close()

	header = "Peak\tLength" + "\t" + "\t".join( sorted( ppeakDict[list(ppeakDict.keys())[0]] ) ) + "\n"
	outwrite = open(outputArg + ".countPerPeak.txt", "w")
	outwrite.write(header)

	for peakKey in sorted(ppeakDict.keys()):
		writeLine = "_".join(peakKey.split("_")[:-1]) + "\t" + peakKey.split("_")[-1]
		writeLine += "\t" + "\t".join( [ str(ppeakDict[peakKey][snpKey]) for snpKey in sorted( list(ppeakDict[peakKey].keys()) ) ] ) + "\n"
		outwrite.write(writeLine)
	outwrite.close()


def motif_centered_peaks(bdtArg, tfPeakArg, offsetArg, testSizeArg):
	openTf = open(tfPeakArg, "r")
	tfLines = openTf.readlines()
	openTf.close()

	tfDict = dict()
	for tfLine in tfLines:
		tfLine = tfLine.rstrip().split()
		tfDict[tfLine[3]] = tfLine

	openOffset = open(offsetArg, "r")
	offsetLines = openOffset.readlines()
	openOffset.close()

	resize = int(int(testSizeArg)/2)

	for offsetLine in offsetLines:
		offsetLine = offsetLine.rstrip().split()
		peakName = offsetLine[3]

		coord_point = int((int(offsetLine[1]) + int(offsetLine[2]))/2) + (int(offsetLine[1]) + int(offsetLine[2]))%2 # for bed inputs
		start = coord_point - 1 - resize
		end = coord_point + resize - 2 #- 1
		tfDict[peakName][1] = str(start)
		tfDict[peakName][2] = str(end)

	outwrite = open(tfPeakArg + ".motifCentered.width" + testSizeArg + "__tmp__", "w")
	for tfKey in sorted(tfDict.keys()):
		writeLine = "\t".join(tfDict[tfKey]) + "\n"
		outwrite.write(writeLine)
	outwrite.close()

	cmd_sort = bdtArg + " sort"
	cmd_sort += " -g " + "~/dnd/genome/GRCh38.p14.genome.chrs.bed"
	cmd_sort += " -i " + tfPeakArg + ".motifCentered.width" + testSizeArg + "__tmp__"
	cmd_sort += " > " + tfPeakArg + ".motifCentered.width" + testSizeArg

	print (cmd_sort, "\n")
	subprocess.getoutput(cmd_sort)

	cmd_rm = remove_temp( "/".join(tfPeakArg.split("/")[:-1]), tfPeakArg.split("/")[-1] + ".motifCentered.width" + testSizeArg + "__tmp__" )
	subprocess.getoutput(cmd_rm)
	#subprocess.getoutput("rm " + tfPeakArg + ".motifCentered.width" + testSizeArg + "__tmp__")


def motifCenteredCounts(peakArg, vcfArg, bkgdPeakArg, bkgdVcfArg, outputArg, sizeArg, sampArg, editArg, editVars):
	range_start = int(0 - int(sizeArg)/2)
	range_end = int(int(sizeArg)/2 +1)

	for r in range(10):
		counts = dict()
		for dis in range(range_start, range_end):
			counts[dis] = [0,0]

		for n in [1,2]:
			if n == 1:
				peakFile = peakArg
				vcfFile = vcfArg
			elif n == 2:
				peakFile = bkgdPeakArg
				vcfFile = bkgdVcfArg

			openPeak = open(peakFile, "r")
			peakLines = openPeak.readlines()
			openPeak.close()

			openVcf = open(vcfFile, "r")
			vcfLines = openVcf.readlines()
			openVcf.close()
			vcfLines = [x for x in vcfLines if x[0] != "#"]

			randIdx = sorted(random.sample(range(len(peakLines)), sampArg))
			peakLines_sampled = [peakLines[i] for i in randIdx]

			newLines = list()
			for vcfLine in vcfLines:
				testLine = vcfLine.rstrip().split()
				vcfChr = testLine[0]
				vcfPos = int(testLine[1])
				vcfRef = testLine[3]
				vcfAlt = testLine[4]
				vcfCount = int(testLine[-1].split(",")[-1])
				vcfDp = int(testLine[-1].split(":")[0])
				vcfFrac = vcfCount/vcfDp ###

				if editArg == "DnD":
					if vcfRef + vcfAlt not in editVars:	continue
				elif editArg == "nonDnD":
					if vcfRef + vcfAlt in editVars:	continue

				for peakLine in peakLines_sampled:
					peakLine = peakLine.rstrip().split()
					peak_chr = peakLine[0]
					peak_start = int(peakLine[1]) +1
					peak_end = int(peakLine[2]) +1
					peak_name = peakLine[3]
					peak_center = peak_start +int(sizeArg)/2 -1

					if peak_chr != vcfChr:	continue
					elif vcfPos < peak_start:	break
					elif vcfPos >= peak_start and vcfPos <= peak_end:
						distance = vcfPos - peak_center
						counts[distance][n-1] += vcfCount

		outwrite = open(outputArg + "_sampled_n" + str(sampArg) + "_" + str(r+1) + ".txt", "w")
		outwrite.write("Distance\tTF\tBackground\n")
		for disKey in sorted(counts.keys()):
			writeLine = str(disKey) + "\t" + str(counts[disKey][0]) + "\t" + str(counts[disKey][1]) + "\n"
			outwrite.write(writeLine)
		outwrite.close()


def getNt(fastaFilesArg, outputArg, idxArg, fwdArg, bwdArg):
	seqLines = list()
	for fastaFile in fastaFilesArg:
		openFile = open(fastaFile, "r")
		fileLines = openFile.readlines()
		openFile.close()
		fileLines = filter(lambda x: x[0] != ">", fileLines)
		seqLines += fileLines

	countDict = dict()
	for seqLine in seqLines:
		seqLine = seqLine.strip()

		start = idxArg-fwdArg
		if bwdArg:	end = idxArg+bwdArg+1
		else:	end = idxArg+1
		ctx = seqLine[start:end]

		if not ctx in countDict:	countDict[ctx] = 0
		countDict[ctx] += 1

	writeLines = ["Context\tCount\n"]
	writeLines += [x + "\t" + str(countDict[x]) + "\n" for x in sorted(countDict.keys())]

	outwrite = open(outputArg, "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


def revComp(fastaArg):
	openFasta = open(fastaArg, "r")
	fastaLines = openFasta.readlines()
	openFasta.close()

	compDict = {"A":"T", "T":"A", "G":"C", "C":"G", "N":"N"}

	writeLines = list()
	for fastaLine in fastaLines:
		if fastaLine.startswith(">"):
			writeLines.append(fastaLine)
		else:
			seq = fastaLine.strip()
			compSeq = [compDict[nt] for nt in list(seq)]
			compSeq = compSeq[::-1]
			compSeq = "".join(compSeq) + "\n"
			writeLines.append(compSeq)

	if fastaArg.endswith(".fa"):
		outwrite = open(fastaArg.replace(".fa", "_revComp.fa"), "w")
	elif fastaArg.endswith(".fasta"):
		outwrite = open(fastaArg.replace(".fasta", "_revComp.fasta"), "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


if __name__ == '__main__':
	parser = argparse.ArgumentParser('--mode only supports homer2 or sea')
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



