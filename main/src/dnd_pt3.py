
import argparse
import os
import sys
import subprocess
import gzip
import random


###############################
### Step 6. Peak statistics ###
###############################

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

	if not outputArg:
		outputArg = vcfArg[:-3] + "redunFlt.vcf"

	outwrite = open(outputArg, "w")
	redunCheck = list()
	for vcfLine in vcfLines:
		tmp = vcfLine.split()
		tmp_check = "_".join(tmp[0:2] + tmp[3:5])
		if tmp_check not in redunCheck:
			redunCheck.append(tmp_check)
			outwrite.write(vcfLine)
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
			count_vcf = int(vcfLine[-1].split(",")[-1]) ###

			if chrm != chrm_vcf:    continue
			elif start > coord_vcf: continue
			elif end < coord_vcf:   continue
	
			snpKey = ref_vcf + ">" + alt_vcf
			snpKey = revcompDict[snpKey]
			ppeakDict[peakKey][snpKey] += count_vcf ###

		for tmpKey in ppeakDict[peakKey].keys():
			if not tmpKey in countDict:	countDict[tmpKey] = 0.0
			if ppeakDict[peakKey][tmpKey] == 0:	continue
			countDict[tmpKey] += float(ppeakDict[peakKey][tmpKey])/peakrange*100/len(peakLines)

	header = "Ref\tAlt\tSNVs_100bp_peak\n"
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
				vcfCount = int(testLine[-1].split(",")[-1]) ###
				vcfDp = int(testLine[-1].split(":")[0]) ###
				vcfFrac = vcfCount/vcfDp ###

				if editArg == "DnD":
					if vcfRef + vcfAlt not in editVars:
						continue
				elif editArg == "nonDnD":
					if vcfRef + vcfAlt in editVars:
						continue

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
						counts[distance][n-1] += vcfCount ###

		outwrite = open(outputArg + "_sampled_n" + str(sampArg) + "_" + str(r+1) + ".txt", "w")
		outwrite.write("Distance\tTF\tBackground\n")
		for disKey in sorted(counts.keys()):
			writeLine = str(disKey) + "\t" + str(counts[disKey][0]) + "\t" + str(counts[disKey][1]) + "\n"
			outwrite.write(writeLine)
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


def sep_tf_peaks_sea(bdtArg, sampleArg, motifArg, sizeArg, dirArg, outputArg):
	outputArg += "step6_tfpeaks_sea/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	if not os.path.exists(outputArg + sampleArg):
		os.mkdir(outputArg + sampleArg)

	seqFile = dirArg + "step5_peaks/merged/meme_sea/sequences.tsv" #
	siteFile = dirArg + "step5_peaks/merged/meme_sea/sites.tsv" #
	npFile = dirArg + "step5_peaks/merged/peaks_bkflt.narrowPeak" #
	
	openSeq = open(seqFile, "r")
	seqLines = openSeq.readlines()
	openSeq.close()

	seqLines = [x.rstrip().split() for x in seqLines[1:] if x[0] != "#"]
	usepeaks = list()
	for motif in motifArg:
		usepeaks += [x[2] for x in seqLines if x[1] == motif and x[4] == "tp"]
	print (len(usepeaks))

	openSite = open(siteFile, "r")
	siteLines = openSite.readlines()
	openSite.close()

	motifLines = list()
	for motif in motifArg:
		tmpLines = [x for x in siteLines[1:] if x[0] != "#" and x.strip().split()[0] == motif]
		tmpLines = [x for x in tmpLines if x.rstrip().split()[1] in usepeaks]
		motifLines += tmpLines

	siteDict = dict()
	for motifLine in motifLines:
		motifLine = motifLine.rstrip().split()
		if motifLine[1] not in siteDict:
			siteDict[motifLine[1]] = [int(motifLine[2]), int(motifLine[3]), motifLine[4], float(motifLine[5])]
		else:
			if siteDict[motifLine[1]][3] < float(motifLine[5]):
				siteDict[motifLine[1]] = [int(motifLine[2]), int(motifLine[3]), motifLine[4], float(motifLine[5])]
			else:	pass
	usepeaks = list(siteDict.keys())
	
	openNp = open(npFile, "r")
	npLines = openNp.readlines()
	openNp.close()

	tf_npLines = [x for x in npLines if x.rstrip().split()[3] in usepeaks]
	bkgd_npLines = [x for x in npLines if x.rstrip().split()[3] not in usepeaks]

	tf_npDict = dict()
	for tf_npLine in tf_npLines:
		tf_npLine = tf_npLine.rstrip().split()
		tf_npDict[tf_npLine[3]] = [tf_npLine[0], int(tf_npLine[1]), int(tf_npLine[1])] + tf_npLine[3:] # chr start start ...
		center_pos = int( sum( siteDict[tf_npLine[3]][:2] )/2 + sum( siteDict[tf_npLine[3]][:2] )%2 )
		tf_npDict[tf_npLine[3]][1] += center_pos - int(int(sizeArg)/2)
		tf_npDict[tf_npLine[3]][2] += center_pos + int(int(sizeArg)/2) -1 ###
		tf_npDict[tf_npLine[3]][1] = str(tf_npDict[tf_npLine[3]][1])
		tf_npDict[tf_npLine[3]][2] = str(tf_npDict[tf_npLine[3]][2])

	outwrite = open(outputArg + sampleArg + "/motif_TF.narrowPeak", "w")
	outwrite.write("".join(tf_npLines))
	outwrite.close()

	outwrite = open(outputArg + sampleArg + "/motif_bkgd.narrowPeak", "w")
	outwrite.write("".join(bkgd_npLines))
	outwrite.close()

	peakorder = [x.rstrip().split()[3] for x in tf_npLines]
	outwrite = open(outputArg + sampleArg + "/motif_TF.narrowPeak.motifCentered.width" + sizeArg, "w")
	for peak in peakorder:
		outwrite.write("\t".join(tf_npDict[peak]) + "\n")
	outwrite.close()


def sep_tf_peaks_motif(bdtArg, hm2Arg, hm2RefArg, sampleArg, motifArg, motifStgArg, dirArg, outputArg):
	outputArg += "step6_tfpeaks_motif/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	if not os.path.exists(outputArg + sampleArg):
		os.mkdir(outputArg + sampleArg)

	smfile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
	smfile = [dirArg + "step5_peaks/" + sampleArg + "/" + x for x in smfile if x.endswith(".bed")][0]

	motif_peaks = list()
	for i in range(len(motifArg)):
		cmdhm2_find = hm2Arg + " " + smfile
		cmdhm2_find += " " + hm2RefArg
		cmdhm2_find += " " + outputArg + sampleArg + "/homer2_denovo"
		cmdhm2_find += " -find " + motifArg[i]
		cmdhm2_find += " > " +  outputArg + sampleArg + "/target_" + str(i+1) + ".txt"
		motif_peaks.append(outputArg + sampleArg + "/target_" + str(i+1) + ".txt")

		print (cmdhm2_find)
		subprocess.getoutput(cmdhm2_find)
	if motifStgArg:
		for i in range(len(motifStgArg)):
			cmdhm2_find = hm2Arg + " " + smfile_b
			cmdhm2_find += " " + hm2RefArg
			cmdhm2_find += " " + outputArg + sampleArg + "/homer2_denovo"
			cmdhm2_find += " -find " + motifStgArg[i]
			cmdhm2_find += " > " +  outputArg + sampleArg + "/target_bkgd_" + str(i+1) + ".txt"
			motif_peaks.append(outputArg + sampleArg + "/target_bkgd_" + str(i+1) + ".txt")
			
			print (cmdhm2_find)
			subprocess.getoutput(cmdhm2_find)

	if len(motif_peaks) != 1:
		print ("merging de novo motif peaks\n")
		merge_beds(motif_peaks)
	else:
		cmd_mv = "mv " + outputArg + sampleArg + "/target_1.txt " + outputArg + sampleArg + "/target.txt"
		print (cmd_mv, "\n")
		subprocess.getoutput(cmd_mv)

	print ("separating TF and background peaks\n")
	openFile = open(outputArg + sampleArg + "/target.txt", "r")
	fileLines = openFile.readlines()
	openFile.close()

	peakDict = dict()
	for fileLine in fileLines[1:]:
        	fileLine = fileLine.rstrip().split()
	        if fileLine[0] not in peakDict:
        	        peakDict[fileLine[0]] = [float(fileLine[-1]), int(fileLine[1]), fileLine[-2]]
	        else:
        	        if float(fileLine[-1]) > peakDict[fileLine[0]][0]:
                	        peakDict[fileLine[0]] = [float(fileLine[-1]), int(fileLine[1]), fileLine[-2]]

	openSummit = open(smfile, "r")
	summitLines = openSummit.readlines()
	openSummit.close()

	summitDict = dict()
	for summitLine in summitLines:
		summitLine = summitLine.rstrip().split()
		summitDict[summitLine[3]] = [summitLine[0], int(summitLine[1]), int(summitLine[2])] ###

	writeLines = list()
	for peak in sorted(list(peakDict.keys())):
		writeLine = summitDict[peak][0] + "\t" + str(summitDict[peak][1] + peakDict[peak][1]) + "\t" + str(summitDict[peak][2] + peakDict[peak][1]) ###
		writeLine += "\t" + peak + "\t" + str(peakDict[peak][1]) + "\t" + peakDict[peak][2] + "\n" ###
		writeLines.append(writeLine)

	outwrite = open(outputArg + sampleArg + "/__tmp__motif_TF.offset.bed", "w")
	outwrite.write("".join(writeLines))
	outwrite.close()

	cmd_sort = bdtArg + " sort"
	cmd_sort += " -g " + "dnd_acc/bedtools_genome.txt"
	cmd_sort += " -i " + outputArg + sampleArg + "/__tmp__motif_TF.offset.bed"
	cmd_sort += " > " + outputArg + sampleArg + "/motif_TF.offset.bed"
	print (cmd_sort, "\n")
	subprocess.getoutput(cmd_sort)
	subprocess.getoutput("rm " + outputArg + sampleArg + "/__tmp__motif_TF.offset.bed")

	cmd_plt = "app_src/plot_tf_offset.R"
	cmd_plt += " " + outputArg + sampleArg + "/motif_TF.offset.bed"
	subprocess.getoutput(cmd_plt)


	npfile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
	npfile = [x for x in npfile if x.endswith(".narrowPeak")]
	if len(npfile) == 0:
		npfile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
		npfile = [x for x in npfile if x.endswith(".bed")][0]
	else:
		npfile = npfile[0]
	openTarget = open(dirArg + "step5_peaks/" + sampleArg + "/" + npfile, "r")
	targetLines = openTarget.readlines()
	openTarget.close()

	writeTargetLines = list()
	othersLines = list()
	for targetLine in targetLines:
        	if targetLine.split()[3] in peakDict:
                	writeTargetLines.append(targetLine)
	        else:
        	        targetLine = targetLine.rstrip().split()
                	targetLine[3] = targetLine[3] + "_t"
	                targetLine = "\t".join(targetLine) + "\n"
        	        othersLines.append(targetLine)

	outwrite = open(outputArg + sampleArg + "/motif_TF.narrowPeak", "w")
	outwrite.write("".join(writeTargetLines))
	outwrite.close()

	outwrite = open(outputArg + sampleArg + "/motif_bkgd.narrowPeak", "w")
	outwrite.write("".join(othersLines))
	outwrite.close()
	

def non_Dnd_vcf(vcfArg):
	openFile = open(vcfArg, "r")
	fileLines = openFile.readlines()
	openFile.close()

	writeLines = list()
	for fileLine in fileLines:
		if fileLine.startswith("#"):
			writeLines.append(fileLine)
		elif fileLine.split()[3] + fileLine.split()[4] not in ["CT", "GA"]:
			writeLines.append(fileLine)
	
	outwrite = open(vcfArg.replace("vcf", "nonDnD.vcf"), "w")
	outwrite.write("".join(writeLines))
	outwrite.close()


def flt_snvs(vcfArg, varArgs):
	varArgs = varArgs.strip('"')
	varArgs = varArgs.split(",")
	varArgs = [x.strip(" ").replace(">", "") for x in varArgs]

	outVcfArg = vcfArg.replace("SNVs.vcf", "_".join(varArgs) + ".SNVs.vcf")

	snvDict = dict()
	for varArg in varArgs:
		if varArg[0] not in snvDict:
			snvDict[ varArg[0] ] = [ varArg[1] ]
		else:
			snvDict[ varArg[0] ].append(varArg[1])

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
	cmd_sort += " -g " + "dnd_acc/bedtools_genome.txt"
	cmd_sort += " -i " + tfPeakArg + ".motifCentered.width" + testSizeArg + "__tmp__"
	cmd_sort += " > " + tfPeakArg + ".motifCentered.width" + testSizeArg

	print (cmd_sort, "\n")
	subprocess.getoutput(cmd_sort)

	subprocess.getoutput("rm " + tfPeakArg + ".motifCentered.width" + testSizeArg + "__tmp__")


def stats_peaks_sea(bdtArg, v2bArg, sampleArg, testSizeArg, peakSampleSizeArg, dirArg, outputArg, varArgs):
	outputArg += "step6_tfpeaks_sea/" + sampleArg + "/"

	varArgs = varArgs.strip('"')
	varArgs = varArgs.split(",")
	varArgs = [x.strip(" ").replace(">", "") for x in varArgs]

	varSuffix = "_".join(varArgs)
	vcfArg =  dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith("flt.SNVs.vcf")][0]
	try:
		dndVcfArg = dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and varSuffix in x][0]
	except:
		dndVcfArg = flt_snvs(vcfArg, varArgs)
		

	peakFiles = [ outputArg + "motif_TF.narrowPeak", outputArg + "motif_bkgd.narrowPeak", outputArg + "motif_TF.narrowPeak" + ".motifCentered.width" + testSizeArg]
	for peakFile in peakFiles:
		cmdbt_vcf = bdtArg + " intersect -header -wa -a " + dndVcfArg
		cmdbt_vcf += " -b " + peakFile
		cmdbt_vcf += " > " + peakFile + "." + varSuffix + ".vcf"
		print (cmdbt_vcf)
		subprocess.getoutput(cmdbt_vcf)

		cmdbt_vcf = bdtArg + " intersect -header -wa -a " + vcfArg
		cmdbt_vcf += " -b " + peakFile
		cmdbt_vcf += " > " + peakFile + ".vcf"
		print (cmdbt_vcf)
		subprocess.getoutput(cmdbt_vcf)
	
		fltRedunVcfs(peakFile + ".vcf", peakFile + ".vcf")
		non_Dnd_vcf(peakFile + ".vcf")	

		cmdbt_edpk = bdtArg + " intersect -header -wa -a " + peakFile
		cmdbt_edpk += " -b " + peakFile + ".vcf"
		cmdbt_edpk += " > " + peakFile + ".editpeaks"
		print (cmdbt_edpk)
		subprocess.getoutput(cmdbt_edpk)

		cmdbt_uqpk = bdtArg + " merge -i " + peakFile + ".editpeaks"
		cmdbt_uqpk += " -c 4,5,6,7,8,9,10 -o first,first,first,first,first,first,first"
		cmdbt_uqpk += " > " + peakFile
		print (cmdbt_uqpk)
		subprocess.getoutput(cmdbt_uqpk)

		subprocess.getoutput("rm " + peakFile + ".editpeaks")

	for peakFile in peakFiles:
		if not ".width" + testSizeArg in peakFile:
			cmd_resize = "app_src/app_resize_narrowPeak.R"
			cmd_resize += " " + peakFile
			cmd_resize += " " + testSizeArg
			print (cmd_resize)
			resizeRes = subprocess.getoutput(cmd_resize)
			if resizeRes.endswith("Execution halted"):
				print ("Resizing Error; re-running with bed version")
				cmd_resize = "app_src/app_resize.R"
				cmd_resize += " " + peakFile
				cmd_resize += " " + testSizeArg
				cmd_resize += " " + peakFile + ".width" + testSizeArg
				print (cmd_resize)
				resizeRes = subprocess.getoutput(cmd_resize)
			cmdbt_vcf_wd = bdtArg + " intersect -header -wa -a " + peakFile + "." + varSuffix + ".vcf"
			cmdbt_vcf_wd += " -b " + peakFile + ".width" + testSizeArg
			cmdbt_vcf_wd += " > " + peakFile + ".width" + testSizeArg + "." + varSuffix + ".vcf"
			print (cmdbt_vcf_wd)
			subprocess.getoutput(cmdbt_vcf_wd)

			cmdbt_vcf_wd = bdtArg + " intersect -header -wa -a " + peakFile + ".nonDnD.vcf"
			cmdbt_vcf_wd += " -b " + peakFile + ".width" + testSizeArg
			cmdbt_vcf_wd += " > " + peakFile + ".width" + testSizeArg + ".nonDnD.vcf"
			print (cmdbt_vcf_wd)
			subprocess.getoutput(cmdbt_vcf_wd)
			print ("\n")

	if not os.path.exists(outputArg + "summitCenteredCounts"):
		os.mkdir(outputArg + "summitCenteredCounts")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf", #".CT_GA.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf", #".CT_GA.vcf",
				outputArg + "summitCenteredCounts/summitCentered", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "summitCenteredCounts"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)
	
	if not os.path.exists(outputArg + "summitCenteredCounts_nonDnD"):
		os.mkdir(outputArg + "summitCenteredCounts_nonDnD")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "summitCenteredCounts_nonDnD/summitCentered", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "summitCenteredCounts_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt) 

	countOnlySnvs(outputArg + "motif_TF.narrowPeak.vcf", outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "target")
	countOnlySnvs(outputArg + "motif_bkgd.narrowPeak.vcf", outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "others")

	cmdsnr = "app_src/app_snvContext_v2.R"
	cmdsnr += " " + outputArg + "target.countNormalized.txt"
	cmdsnr += " " + outputArg + "others.countNormalized.txt"
	cmdsnr += " " + outputArg
	print (cmdsnr, "\n")
	subprocess.getoutput(cmdsnr)

	if not os.path.exists(outputArg + "motifCenteredCounts"):
		os.mkdir(outputArg + "motifCenteredCounts")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg + "." + varSuffix + ".vcf", #".CT_GA.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf", #".CT_GA.vcf",
				outputArg + "motifCenteredCounts/motifCentered", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "motifCenteredCounts"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + "motifCenteredCounts_nonDnD"):
		os.mkdir(outputArg + "motifCenteredCounts_nonDnD")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motifCenteredCounts_nonDnD/motifCentered", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "motifCenteredCounts_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)


def stats_peaks(bdtArg, v2bArg, sampleArg, testSizeArg, peakSampleSizeArg, dirArg, outputArg, varArgs):
	outputArg += "step6_tfpeaks_motif/" + sampleArg + "/"

	varArgs = varArgs.strip('"')
	varArgs = varArgs.split(",")
	varArgs = [x.strip(" ").replace(">", "") for x in varArgs]
	varSuffix = "_".join(varArgs)

	motif_centered_peaks(bdtArg, outputArg + "motif_TF.narrowPeak", outputArg + "motif_TF.offset.bed", testSizeArg)

	vcfArg =  dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and varSuffix not in x][0]
	dndVcfArg = dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and varSuffix in x][0]

	peakFiles = [ outputArg + "motif_TF.narrowPeak", outputArg + "motif_bkgd.narrowPeak", outputArg + "motif_TF.narrowPeak" + ".motifCentered.width" + testSizeArg]
	for peakFile in peakFiles:
		cmdbt_vcf = bdtArg + " intersect -header -wa -a " + dndVcfArg
		cmdbt_vcf += " -b " + peakFile
		cmdbt_vcf += " > " + peakFile + "." + varSuffix + ".vcf"
		print (cmdbt_vcf)
		subprocess.getoutput(cmdbt_vcf)
	
		cmdbt_vcf = bdtArg + " intersect -header -wa -a " + vcfArg
		cmdbt_vcf += " -b " + peakFile
		cmdbt_vcf += " > " + peakFile + ".vcf"
		print (cmdbt_vcf)
		subprocess.getoutput(cmdbt_vcf)
	
		fltRedunVcfs(peakFile + ".vcf", peakFile + ".vcf")

		non_Dnd_vcf(peakFile + ".vcf")

		cmdbt_edpk = bdtArg + " intersect -header -wa -a " + peakFile
		cmdbt_edpk += " -b " + peakFile + ".vcf"
		cmdbt_edpk += " > " + peakFile + ".editpeaks"
		print (cmdbt_edpk)
		subprocess.getoutput(cmdbt_edpk)
		
		cmdbt_uqpk = bdtArg + " merge -i " + peakFile + ".editpeaks"
		cmdbt_uqpk += " -c 4,5,6,7,8,9,10 -o first,first,first,first,first,first,first"
		cmdbt_uqpk += " > " + peakFile
		print (cmdbt_uqpk)
		subprocess.getoutput(cmdbt_uqpk)

		subprocess.getoutput("rm " + peakFile + ".editpeaks")
		
	
	for peakFile in peakFiles: ### 2024.07.09 end
		if not ".width" + testSizeArg in peakFile:
			cmd_resize = "app_src/app_resize_narrowPeak.R"
			cmd_resize += " " + peakFile
			cmd_resize += " " + testSizeArg
			print (cmd_resize)
			resizeRes = subprocess.getoutput(cmd_resize)
			if resizeRes.endswith("Execution halted"):
				print ("Resizing Error; re-running with bed version")
				cmd_resize = "app_src/app_resize.R"
				cmd_resize += " " + peakFile
				cmd_resize += " " + testSizeArg
				cmd_resize += " " + peakFile + ".width" + testSizeArg
				print (cmd_resize)
				resizeRes = subprocess.getoutput(cmd_resize)

			cmdbt_vcf_wd = bdtArg + " intersect -header -wa -a " + peakFile + "." + varSuffix + ".vcf"
			cmdbt_vcf_wd += " -b " + peakFile + ".width" + testSizeArg
			cmdbt_vcf_wd += " > " + peakFile + ".width" + testSizeArg + "." + varSuffix + ".vcf"
			print (cmdbt_vcf_wd)
			subprocess.getoutput(cmdbt_vcf_wd)

			cmdbt_vcf_wd = bdtArg + " intersect -header -wa -a " + peakFile + ".nonDnD.vcf"
			cmdbt_vcf_wd += " -b " + peakFile + ".width" + testSizeArg
			cmdbt_vcf_wd += " > " + peakFile + ".width" + testSizeArg + ".nonDnD.vcf"
			print (cmdbt_vcf_wd)
			subprocess.getoutput(cmdbt_vcf_wd)
			print ("\n")

	if not os.path.exists(outputArg + "summitCenteredCounts"):
		os.mkdir(outputArg + "summitCenteredCounts")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf", 
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf",
				outputArg + "summitCenteredCounts/summitCentered", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "summitCenteredCounts"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + "summitCenteredCounts_nonDnD"):
		os.mkdir(outputArg + "summitCenteredCounts_nonDnD")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "summitCenteredCounts_nonDnD/summitCentered", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)	
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "summitCenteredCounts_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	countOnlySnvs(outputArg + "motif_TF.narrowPeak.vcf", outputArg + "motif_TF.narrowPeak.width" + testSizeArg, outputArg + "target")
	countOnlySnvs(outputArg + "motif_bkgd.narrowPeak.vcf", outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "others")

	cmdsnr = "app_src/app_snvContext_v2.R" 
	cmdsnr += " " + outputArg + "target.countNormalized.txt"
	cmdsnr += " " + outputArg + "others.countNormalized.txt"
	cmdsnr += " " + outputArg
	print (cmdsnr, "\n")
	subprocess.getoutput(cmdsnr)


	if not os.path.exists(outputArg + "motifCenteredCounts"):
		os.mkdir(outputArg + "motifCenteredCounts")	
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg + "." + varSuffix + ".vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + "." + varSuffix + ".vcf",
				outputArg + "motifCenteredCounts/motifCentered", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "motifCenteredCounts"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + "motifCenteredCounts_nonDnD"):
		os.mkdir(outputArg + "motifCenteredCounts_nonDnD")
	motifCenteredCounts(outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg, outputArg + "motif_TF.narrowPeak.motifCentered.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + "motif_bkgd.narrowPeak.width" + testSizeArg + ".nonDnD.vcf",
				outputArg + "motifCenteredCounts_nonDnD/motifCentered", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + "motifCenteredCounts_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)


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
                if bwdArg:      end = idxArg+bwdArg+1
                else:   end = idxArg+1
                ctx = seqLine[start:end]

                if not ctx in countDict:
                        countDict[ctx] = 0
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

def sep_tf_peaks_chip(bdtArg, v2bArg, targetPeakArg, sampleArg, testSizeArg, peakSampleSizeArg, dirArg, outputArg, varArgs):
	outputArg += "step6_tfpeaks_chipseq/"
	if not os.path.exists(outputArg):
		os.mkdir(outputArg)

	if not os.path.exists(outputArg + sampleArg):
		os.mkdir(outputArg + sampleArg)

	varArgs = varArgs.strip('"')
	varArgs = varArgs.split(",")
	varArgs = [x.strip(" ").replace(">", "") for x in varArgs]


	npfile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
	npfile = [x for x in npfile if x.endswith(".narrowPeak")][0]
	cmdbt_inter = bdtArg + " intersect -wa"
	cmdbt_inter += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + npfile
	cmdbt_inter += " -b " + targetPeakArg
	cmdbt_inter += " | " + bdtArg + " sort -i -"
	cmdbt_inter += " > " + outputArg + sampleArg + "/target.narrowPeak"
	print (cmdbt_inter)
	subprocess.getoutput(cmdbt_inter)

	print ("fltRedunPeaks")
	fltRedunPeaks(outputArg + sampleArg + "/target.narrowPeak")
	
	smfile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
	smfile = [x for x in smfile if x.endswith(".bed")][0]
	cmdbt_intersm =  bdtArg + " intersect -wa"
	cmdbt_intersm += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + smfile
	cmdbt_intersm += " -b " + outputArg + sampleArg + "/target.narrowPeak"
	cmdbt_intersm += " > " + outputArg + sampleArg + "/target.summits.bed"
	print (cmdbt_intersm)
	subprocess.getoutput(cmdbt_intersm)

	cmdresize = "app_src/app_resize_narrowPeak.R"
	cmdresize += " " + outputArg + sampleArg + "/target.narrowPeak"
	cmdresize += " " + testSizeArg
	print (cmdresize)
	subprocess.getoutput(cmdresize)

	vcfFile = [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and "CT_GA" not in x][0]
	cmdbt_vcf = bdtArg + " intersect -wa -a " + dirArg + "step4_snvs/" + sampleArg + "/" + vcfFile
	cmdbt_vcf += " -b " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg
	cmdbt_vcf += " > " + outputArg + sampleArg + "/__tmp__target.narrowPeak.width" + testSizeArg + ".vcf"
	print (cmdbt_vcf)
	subprocess.getoutput(cmdbt_vcf)

	fltRedunVcfs(outputArg + sampleArg + "/__tmp__target.narrowPeak.width" + testSizeArg + ".vcf", outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf")

	countOnlySnvs(outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf", outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/target")
	print ("\n")

	cmdbt_inter = bdtArg + " subtract -A"
	cmdbt_inter += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + npfile
	cmdbt_inter += " -b " + targetPeakArg
	cmdbt_inter += " | " + bdtArg + " subtract -A -a - -b " + outputArg + sampleArg + "/target.narrowPeak"
	cmdbt_inter += " | " + bdtArg + " sort -i - > " + outputArg + sampleArg + "/others.narrowPeak"
	print (cmdbt_inter)
	subprocess.getoutput(cmdbt_inter)

	fltRedunPeaks(outputArg + sampleArg + "/others.narrowPeak")

	cmdbt_intersm =  bdtArg + " intersect -wa"
	cmdbt_intersm += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + smfile
	cmdbt_intersm += " -b " + outputArg + sampleArg + "/others.narrowPeak"
	cmdbt_intersm += " > " + outputArg + sampleArg + "/others.summits.bed"
	print (cmdbt_intersm)
	subprocess.getoutput(cmdbt_intersm)

	cmdresize = "app_src/app_resize_narrowPeak.R"
	cmdresize += " " + outputArg + sampleArg + "/others.narrowPeak"
	cmdresize += " " + testSizeArg
	print (cmdresize)
	subprocess.getoutput(cmdresize)       

	cmdbt_vcf = bdtArg + " intersect -wa -a " + dirArg + "step4_snvs/" + sampleArg + "/" + vcfFile
	cmdbt_vcf += " -b " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg
	cmdbt_vcf += " > " + outputArg + sampleArg + "/__tmp__others.narrowPeak.width" + testSizeArg + ".vcf"
	print (cmdbt_vcf)
	subprocess.getoutput(cmdbt_vcf)

	fltRedunVcfs(outputArg + sampleArg + "/__tmp__others.narrowPeak.width" + testSizeArg + ".vcf", outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf")
	
	countOnlySnvs(outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf", outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/others")
	cmd_rm = "rm " + outputArg + sampleArg + "/__tmp__*"
	subprocess.getoutput(cmd_rm)

	cmdsnr = "app_src/app_snvContext_v2.R " + outputArg + sampleArg + "/target.countNormalized.txt"
	cmdsnr += " " + outputArg + sampleArg + "/others.countNormalized.txt"
	cmdsnr += " " + outputArg + sampleArg
	print (cmdsnr, "\n")
	subprocess.getoutput(cmdsnr)

	cmdv2b_tf = v2bArg + " < " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf" + " > " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.bed"
	print (cmdv2b_tf)
	subprocess.getoutput(cmdv2b_tf)

	cmdv2b_ot = v2bArg + " < " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf" + " > " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.bed"
	print (cmdv2b_ot)
	subprocess.getoutput(cmdv2b_ot)

	### motify to analyze other variants ###
	cmd_wd_tf = "app_src/app_extend_xbp.R " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.bed" + " " + str(int(int(testSizeArg)/2))
	subprocess.getoutput(cmd_wd_tf)
	cmd_wd_ot = "app_src/app_extend_xbp.R " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.bed" + " " + str(int(int(testSizeArg)/2))
	subprocess.getoutput(cmd_wd_ot)

	cmd_bdt_tf = bdtArg + " getfasta -fi dnd_acc/genome.fa"
	cmd_bdt_tf += " -fo " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa"
	cmd_bdt_tf += " -bed " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.bed"
	print (cmd_bdt_tf)
	subprocess.getoutput(cmd_bdt_tf)

	cmd_bdt_tf = bdtArg + " getfasta -fi dnd_acc/genome.fa"
	cmd_bdt_tf += " -fo " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.fa"
	cmd_bdt_tf += " -bed " + outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.bed"
	print (cmd_bdt_tf)
	subprocess.getoutput(cmd_bdt_tf)

	cmd_bdt_ot = bdtArg + " getfasta -fi dnd_acc/genome.fa"
	cmd_bdt_ot += " -fo " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa"
	cmd_bdt_ot += " -bed " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.bed"
	print (cmd_bdt_ot)
	subprocess.getoutput(cmd_bdt_ot)

	cmd_bdt_ot = bdtArg + " getfasta -fi dnd_acc/genome.fa"
	cmd_bdt_ot += " -fo " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.fa"
	cmd_bdt_ot += " -bed " + outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.bed"
	print (cmd_bdt_ot, "\n")
	subprocess.getoutput(cmd_bdt_ot)

	revComp(outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.fa")
	revComp(outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp.fa")


	if not os.path.exists(outputArg + sampleArg + "/summitCenteredCounts"):
		os.mkdir(outputArg + sampleArg + "/summitCenteredCounts")
	motifCenteredCounts(outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf",
				outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf",
				outputArg + sampleArg + "/summitCenteredCounts/summitCentered", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + sampleArg + "/summitCenteredCounts"
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + sampleArg + "/summitCenteredCounts_nonDnD"):
		os.mkdir(outputArg + sampleArg + "/summitCenteredCounts_nonDnD")
	motifCenteredCounts(outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf",
				outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf",
				outputArg + sampleArg + "/summitCenteredCounts_nonDnD/summitCentered", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "app_src/plot_motifCenteredCounts_encode.R " + outputArg + sampleArg + "/summitCenteredCounts_nonDnD"
	subprocess.getoutput(cmd_plt)


	if not os.path.exists(outputArg + sampleArg + "/context_di"):
		os.mkdir(outputArg + sampleArg + "/context_di")
	getNt([outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa", 
		outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp_revComp.fa",],
		outputArg + sampleArg + "/context_di/target.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd1.txt", 100, 1, 0)
	getNt([outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa", 
		outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp_revComp.fa"],
		outputArg + sampleArg + "/context_di/others.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd1.txt", 100, 1, 0)

	if not os.path.exists(outputArg + sampleArg + "/context_tri"):
		os.mkdir(outputArg + sampleArg + "/context_tri")
	getNt([outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa", 
		outputArg + sampleArg + "/target.narrowPeak.width" + testSizeArg + ".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp_revComp.fa",],
		outputArg + sampleArg + "/context_tri/target.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd2.txt", 100, 2, 0)
	getNt([outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg + ".vcf.CtoT_+-" + str(int(int(testSizeArg)/2)) + "bp.fa", 
		outputArg + sampleArg + "/others.narrowPeak.width" + testSizeArg +".vcf.GtoA_+-" + str(int(int(testSizeArg)/2)) + "bp_revComp.fa"],
		outputArg + sampleArg + "/context_tri/others.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd2.txt", 100, 2, 0)

	cmd_plt = "app_src/plot_context.R " + outputArg + sampleArg + "/context_di/target.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd1.txt 2"
	subprocess.getoutput(cmd_plt)
	cmd_plt = "app_src/plot_context.R " + outputArg + sampleArg + "/context_di/others.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd1.txt 2"
	subprocess.getoutput(cmd_plt)

	cmd_plt = "app_src/plot_context.R " + outputArg + sampleArg + "/context_tri/target.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd2.txt 3"
	subprocess.getoutput(cmd_plt)
	cmd_plt = "app_src/plot_context.R " + outputArg + sampleArg + "/context_tri/others.narrowPeak.width" + testSizeArg + ".vcf.CT_GA_+-" + str(int(int(testSizeArg)/2)) + "bp.fwd2.txt 3"
	subprocess.getoutput(cmd_plt)
###


def main(args):
	bdt = subprocess.getoutput("which bedtools")
	mc2 = subprocess.getoutput("which macs2")
	v2b = subprocess.getoutput("which vcf2bed")
	hm2 = subprocess.getoutput("which findMotifsGenome.pl")

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
		if args.mode == "homer2":
			sep_tf_peaks_motif(bdt, hm2, args.hm2ref, sample, args.motif, False, args.Dir, args.Output)
			stats_peaks(bdt, v2b, sample, args.size, args.rand, args.Dir, args.Output, args.variants)
		elif args.mode == "chip":
			sep_tf_peaks_chip(bdt, v2b, args.chipseq, sample, args.size, args.rand, args.Dir, args.Output, args.variants)
		elif args.mode == "sea":
			sep_tf_peaks_sea(bdt, sample, args.motif, args.size, args.Dir, args.Output)
			stats_peaks_sea(bdt, v2b, sample, args.size, args.rand, args.Dir, args.Output, args.variants)


if __name__ == '__main__':
	parser = argparse.ArgumentParser('--mode only supports homer2 or sea')
	parser.add_argument('-d', '--Dir', help = 'directory path', required = True)
	parser.add_argument('-o', '--Output', help = '[Global] (*optional) output directory path', required = False)

	parser.add_argument('--size', help = '[Global] (*optional) test peak width size; default is 200', type = str, default = "200", required = False)
	parser.add_argument('--rand', help = '[Global] (*optional) down-sample peak number; default is 200', type = int, default = 200, required = False)
	parser.add_argument('--sample', help = '[Global] "sample name" or "all" for all samples in <step5>', required = True)
	parser.add_argument('--var', dest = 'variants', help = '[Global] (*optional) expected D&D variants; default is "C>T,G>A"', default = "C>T,G>A", required = False)

	parser.add_argument('--mode', help = '[Global] which mode: "chip", "homer2" or "sea"', nargs = '?', choices = ['chip', 'homer2', 'sea'], default = 'sea', required = True)

	parser.add_argument('--chipseq', help = '[Step 6, --mode:chip] chip-seq reference', required = False)
	parser.add_argument('--homer-ref', dest = 'hm2ref', help = '[Step 6, --mode:homer2] (*optional) homer2 reference; default is "grch38_crgatac"', default = 'grch38_crgatac', required = False)
	parser.add_argument('--motif', help = '[Step 6, --mode:homer2 or sea] <path to the homer2 motif file> for "homer2" or <TF name> for "sea"', nargs = '+', required = False)

	args = parser.parse_args()
	main(args)


