
import argparse
import os
import sys
import subprocess
import gzip
import random
from main.dnd_acc import *


def tf_peaks(modeArg, bdtArg, hm2Arg, hm2RefArg, sampleArg, motifArg, targetPeakArg, sizeArg, dirArg, outputArg, indivArg):
	if modeArg == "sea":
		outputArg += "step6_tfpeaks_sea/"
	elif modeArg == "homer2":
		outputArg += "step6_tfpeaks_homer2/"
	elif modeArg == "chip":
		outputArg += "step6_tfpeaks_chipseq/"
	if not os.path.exists(outputArg):	os.mkdir(outputArg)
	if not os.path.exists(outputArg + sampleArg):	os.mkdir(outputArg + sampleArg)

	if modeArg == "sea":
		if indivArg:
			seqFile = dirArg + "step5_peaks/" + sampleArg + "/meme_sea/sequences.tsv"
			siteFile = dirArg + "step5_peaks/"+ sampleArg + "/meme_sea/sites.tsv"
		else:
			seqFile = dirArg + "step5_peaks/merged/meme_sea/sequences.tsv"
			siteFile = dirArg + "step5_peaks/merged/meme_sea/sites.tsv"

		openSeq = open(seqFile, "r")
		seqLines = openSeq.readlines()
		openSeq.close()
		
		seqLines = [x.rstrip().split("\t") for x in seqLines[1:] if x[0] != "#"]
		usepeaks = list()
		for motif in motifArg:	usepeaks += [x[3] for x in seqLines if x[1] == motif and x[5] == "tp"]
		print (len(usepeaks))

		openSite = open(siteFile, "r")
		siteLines = openSite.readlines()
		openSite.close()

		motifLines = list()
		for motif in motifArg:
			tmpLines = [x for x in siteLines[1:] if x[0] != "#" and x.strip().split()[0] == motif]
			tmpLines = [x for x in tmpLines if x.rstrip().split("\t")[2] in usepeaks]
			motifLines += tmpLines

		peakDict = dict()
		for motifLine in motifLines:
			motifLine = motifLine.rstrip().split()
			if motifLine[2] not in peakDict:
				# key: peak name, items: site start, site end, strand, site score
				peakDict[motifLine[2]] = [int(motifLine[3]), int(motifLine[4]), motifLine[5], float(motifLine[6])]
			else:
				if peakDict[motifLine[2]][3] < float(motifLine[6]):
					peakDict[motifLine[2]] = [int(motifLine[3]), int(motifLine[4]), motifLine[5], float(motifLine[6])]
				else:	pass
		usepeaks = list(peakDict.keys())

	elif modeArg == "homer2":
		smFile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
		smFile = [dirArg + "step5_peaks/" + sampleArg + "/" + x for x in smFile if x.endswith(".bed")][0]

		motif_peaks = list()
		for i in range(len(motifArg)):
			cmd_hm2_find = hm2Arg + " " + smFile
			cmd_hm2_find += " " + hm2RefArg
			cmd_hm2_find += " " + outputArg + sampleArg + "/homer2_denovo"
			cmd_hm2_find += " -find " + motifArg[i]
			cmd_hm2_find += " > " +  outputArg + sampleArg + "/target_" + str(i+1) + ".txt"
			motif_peaks.append(outputArg + sampleArg + "/target_" + str(i+1) + ".txt")
			print (cmd_hm2_find)
			subprocess.getoutput(cmd_hm2_find)
		if len(motif_peaks) != 1:
			print ("merging de novo motif peaks\n")
			merge_beds(motif_peaks)
		else:
			cmd_mv = "mv " + outputArg + sampleArg + "/target_1.txt " + outputArg + sampleArg + "/target.txt"
			print (cmd_mv, "\n")
			subprocess.getoutput(cmd_mv)

		openFile = open(outputArg + sampleArg + "/target.txt", "r")
		motifLines = openFile.readlines()
		openFile.close()

		peakDict = dict()
		for motifLine in motifLines[1:]:
			motifLine = motifLine.rstrip().split()
			if motifLine[0] not in peakDict:
				# key: peak name, items: motif score, offset, strand
				peakDict[motifLine[0]] = [float(motifLine[-1]), int(motifLine[1]), motifLine[-2]]
			else:
				if peakDict[motifLine[0]][0] < float(motifLine[-1]):
					peakDict[motifLine[0]] = [float(motifLine[-1]), int(motifLine[1]), motifLine[-2]]

		openSummit = open(smFile, "r")
		summitLines = openSummit.readlines()
		openSummit.close()
		summitDict = dict()
		for summitLine in summitLines:
			summitLine = summitLine.rstrip().split()
			summitDict[summitLine[3]] = [summitLine[0], int(summitLine[1]), int(summitLine[2])]

		writeLines = list()
		for peak in sorted(list(peakDict.keys())):
			writeLine = summitDict[peak][0] + "\t" + str(summitDict[peak][1] + peakDict[peak][1]) + "\t" + str(summitDict[peak][2] + peakDict[peak][1]) ###
			writeLine += "\t" + peak + "\t" + str(peakDict[peak][1]) + "\t" + peakDict[peak][2] + "\n"
			writeLines.append(writeLine)

		outwrite = open(outputArg + sampleArg + "/__tmp__motif_TF.offset.bed", "w")
		outwrite.write("".join(writeLines))
		outwrite.close()

		cmd_sort = bdtArg + " sort"
		cmd_sort += " -g " + "~/DnD/genome/GRCh38.p14.genome.chrs.bed"
		cmd_sort += " -i " + outputArg + sampleArg + "/__tmp__motif_TF.offset.bed"
		cmd_sort += " > " + outputArg + sampleArg + "/motif_TF.offset.bed"
		print (cmd_sort, "\n")
		subprocess.getoutput(cmd_sort)
		
		cmd_rm = remove_temp(outputArg + sampleArg, "__tmp__")
		subprocess.getoutput(cmd_rm)

		cmd_plt = "Rscript ~/DnD/src/main/tools.R plot_offset"
		cmd_plt += " " + outputArg + sampleArg + "/motif_TF.offset.bed"
		subprocess.getoutput(cmd_plt)

	elif modeArg == "chip":
		macsFiles = os.listdir(dirArg + "step5_peaks/" + sampleArg)
		npFile = [x for x in macsFiles if x.endswith(".narrowPeak")][0]
		smFile = [x for x in macsFiles if x.endswith(".bed")][0]

		cmd_bt_ta_np = bdtArg + " intersect -wa"
		cmd_bt_ta_np += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + npFile
		cmd_bt_ta_np += " -b " + targetPeakArg
		cmd_bt_ta_np += " | " + bdtArg + " sort -i -"
		cmd_bt_ta_np += " > " + outputArg + sampleArg + "/motif_TF.narrowPeak"
		print (cmd_bt_ta_np)
		subprocess.getoutput(cmd_bt_ta_np)
		fltRedunPeaks(outputArg + sampleArg + "/motif_TF.narrowPeak")

		cmd_bt_ta_sm =  bdtArg + " intersect -wa"
		cmd_bt_ta_sm += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + smFile
		cmd_bt_ta_sm += " -b " + outputArg + sampleArg + "/motif_TF.narrowPeak"
		cmd_bt_ta_sm += " > " + outputArg + sampleArg + "/motif_TF.summits.bed"
		print (cmd_bt_ta_sm)
		subprocess.getoutput(cmd_bt_ta_sm)
		fltRedunPeaks(outputArg + sampleArg + "/motif_TF.summits.bed")

		cmd_resize = "Rscript ~/DnD/src/main/tools.R app_resize"

		cmd_resize_np = cmd_resize + " " + outputArg + sampleArg + "/motif_TF.narrowPeak"
		cmd_resize_np += " " + sizeArg + " np"
		print (cmd_resize_np)
		subprocess.getoutput(cmd_resize_np)
		cmd_resize_sm = cmd_resize + " " + outputArg + sampleArg + "/motif_TF.summits.bed"
		cmd_resize_sm += " " + sizeArg + " sm"
		print (cmd_resize_sm)
		subprocess.getoutput(cmd_resize_sm)

		###
		cmd_bt_bk_np = bdtArg + " subtract -A"
		cmd_bt_bk_np += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + npFile
		cmd_bt_bk_np += " -b " + targetPeakArg
		cmd_bt_bk_np += " | " + bdtArg + " subtract -A -a - -b " + outputArg + sampleArg + "/motif_TF.narrowPeak"
		cmd_bt_bk_np += " | " + bdtArg + " sort -i - > " + outputArg + sampleArg + "/motif_bkgd.narrowPeak"
		print (cmd_bt_bk_np)
		subprocess.getoutput(cmd_bt_bk_np)
		fltRedunPeaks(outputArg + sampleArg + "/motif_bkgd.narrowPeak")

		cmd_bt_bk_sm =  bdtArg + " intersect -wa"
		cmd_bt_bk_sm += " -a " + dirArg + "step5_peaks/" + sampleArg + "/" + smFile
		cmd_bt_bk_sm += " -b " + outputArg + sampleArg + "/motif_bkgd.narrowPeak"
		cmd_bt_bk_sm += " > " + outputArg + sampleArg + "/motif_bkgd.summits.bed"
		print (cmd_bt_bk_sm)
		subprocess.getoutput(cmd_bt_bk_sm)
		fltRedunPeaks(outputArg + sampleArg + "/motif_bkgd.summits.bed")

		cmd_resize_np = cmd_resize + " " + outputArg + sampleArg + "/motif_bkgd.narrowPeak"
		cmd_resize_np += " " + sizeArg + " np"
		print (cmd_resize_np)
		subprocess.getoutput(cmd_resize_np)
		cmd_resize_sm = cmd_resize + " " + outputArg + sampleArg + "/motif_bkgd.summits.bed"
		cmd_resize_sm += " " + sizeArg + " sm"
		print (cmd_resize_sm)
		subprocess.getoutput(cmd_resize_sm)

	if modeArg in ["sea", "homer2"]:
		npFile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
		npFile = [x for x in npFile if x.endswith(".narrowPeak")]
		if len(npFile) == 0:
			npFile = os.listdir(dirArg + "step5_peaks/" + sampleArg)
			npFile = [x for x in npFile if x.endswith(".bed")][0]
		else:
			npFile = npFile[0]
	
		openNp = open(dirArg + "step5_peaks/" + sampleArg + "/" + npFile, "r")
		npLines = openNp.readlines()
		openNp.close()

		usepeaks = list(peakDict.keys())
	
		tf_npLines = [x for x in npLines if x.rstrip().split()[3] in usepeaks]
		bkgd_npLines = [x for x in npLines if x.rstrip().split()[3] not in usepeaks]
		bkgd_npLines = [x.replace("_t", "") for x in bkgd_npLines]

		outwrite = open(outputArg + sampleArg + "/motif_TF.narrowPeak", "w")
		outwrite.write("".join(tf_npLines))
		outwrite.close()

		outwrite = open(outputArg + sampleArg + "/motif_bkgd.narrowPeak", "w")
		outwrite.write("".join(bkgd_npLines))
		outwrite.close()

		if modeArg == "sea":
			tf_npDict = dict()
			for tf_npLine in tf_npLines:
				tf_npLine = tf_npLine.rstrip().split()
				tf_npDict[tf_npLine[3]] = [tf_npLine[0], int(tf_npLine[1]), int(tf_npLine[1])] + tf_npLine[3:]
				center_pos = int( sum( peakDict[tf_npLine[3]][:2] )/2 + sum( peakDict[tf_npLine[3]][:2] )%2 )
				tf_npDict[tf_npLine[3]][1] += center_pos - int(int(sizeArg)/2)
				tf_npDict[tf_npLine[3]][2] += center_pos + int(int(sizeArg)/2) -1
				tf_npDict[tf_npLine[3]][1] = str(tf_npDict[tf_npLine[3]][1])
				tf_npDict[tf_npLine[3]][2] = str(tf_npDict[tf_npLine[3]][2])

			peakorder = [x.rstrip().split()[3] for x in tf_npLines]
			outwrite = open(outputArg + sampleArg + "/motif_TF.narrowPeak.motifCentered.width" + sizeArg, "w")
			for peak in peakorder:
				outwrite.write("\t".join(tf_npDict[peak]) + "\n")
			outwrite.close()
		elif modeArg == "homer2":
			motif_centered_peaks(bdtArg, outputArg + sampleArg + "/motif_TF.narrowPeak", outputArg + sampleArg + "/motif_TF.offset.bed", sizeArg)


def stats_peaks(modeArg, bdtArg, v2bArg, sampleArg, testSizeArg, peakSampleSizeArg, dirArg, outputArg, varArgs):
	if modeArg == "sea":
		outputArg += "step6_tfpeaks_sea/"
	elif modeArg == "homer2":
		outputArg += "step6_tfpeaks_homer2/"
	elif modeArg == "chip":
		outputArg += "step6_tfpeaks_chipseq/"
	if not os.path.exists(outputArg):	os.mkdir(outputArg)
	if not os.path.exists(outputArg + sampleArg):	os.mkdir(outputArg + sampleArg)

	varArgs = varArgs.strip('"')
	varArgs = varArgs.split(",")
	varArgs = [x.strip(" ").replace(">", "") for x in varArgs]
	varSuffix = "_".join(varArgs)

	vcfArg =  dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and varSuffix not in x][0]
	try:	dndVcfArg = dirArg + "step4_snvs/" + sampleArg + "/" + [x for x in os.listdir(dirArg + "step4_snvs/" + sampleArg) if x.endswith(".SNVs.vcf") and varSuffix in x][0]
	except:	dndVcfArg = dnd_vcf(vcfArg, varArgs)

	if modeArg in ["sea", "homer2"]:
		peakFiles = [ outputArg + sampleArg + "/motif_TF.narrowPeak", 
				outputArg + sampleArg + "/motif_TF.narrowPeak" + ".motifCentered.width" + testSizeArg,
				outputArg + sampleArg + "/motif_bkgd.narrowPeak"]
	elif modeArg == "chip":
		peakFiles = [ outputArg + sampleArg + "/motif_TF.narrowPeak.width200",
				outputArg + sampleArg + "/motif_TF.summits.bed.width200",
				outputArg + sampleArg + "/motif_bkgd.narrowPeak.width200",
				outputArg + sampleArg + "/motif_bkgd.summits.bed.width200"]
	for peakFile in peakFiles:
		cmd_bt_vcf = bdtArg + " intersect -header -wa -a " + dndVcfArg
		cmd_bt_vcf += " -b " + peakFile
		cmd_bt_vcf += " > " + peakFile + "." + varSuffix + ".vcf"
		print (cmd_bt_vcf)
		subprocess.getoutput(cmd_bt_vcf)

		cmd_bt_vcf = bdtArg + " intersect -header -wa -a " + vcfArg
		cmd_bt_vcf += " -b " + peakFile
		cmd_bt_vcf += " > " + peakFile + ".vcf"
		print (cmd_bt_vcf)
		subprocess.getoutput(cmd_bt_vcf)
	
		fltRedunVcfs(peakFile + ".vcf", peakFile + ".vcf")
		non_Dnd_vcf(peakFile + ".vcf", varArgs)	

		cmd_bt_edpk = bdtArg + " intersect -header -wa -a " + peakFile
		cmd_bt_edpk += " -b " + peakFile + ".vcf"
		cmd_bt_edpk += " > " + peakFile + ".editpeaks"
		print (cmd_bt_edpk)
		subprocess.getoutput(cmd_bt_edpk)

		cmd_bt_uqpk = bdtArg + " merge -i " + peakFile + ".editpeaks"

		if "summits" in peakFile:
			cmd_bt_uqpk += " -c 4,5 -o first,first"
		else:	cmd_bt_uqpk += " -c 4,5,6,7,8,9,10 -o first,first,first,first,first,first,first"
		cmd_bt_uqpk += " > " + peakFile
		print (cmd_bt_uqpk)
		subprocess.getoutput(cmd_bt_uqpk)

		subprocess.getoutput("rm " + peakFile + ".editpeaks")

	for peakFile in peakFiles:
		if not ".width" + testSizeArg in peakFile:
			cmd_resize = "Rscript ~/DnD/src/main/tools.R app_resize"
			cmd_resize += " " + peakFile
			cmd_resize += " " + testSizeArg
			cmd_resize += " np"
			print (cmd_resize)
			subprocess.getoutput(cmd_resize)

			cmd_bt_vcf = bdtArg + " intersect -header -wa -a " + peakFile + ".vcf"
			cmd_bt_vcf += " -b " + peakFile + ".width" + testSizeArg
			cmd_bt_vcf += " > " + peakFile + ".width" + testSizeArg + ".vcf"
			print (cmd_bt_vcf)
			subprocess.getoutput(cmd_bt_vcf)

			cmd_bt_vcf = bdtArg + " intersect -header -wa -a " + peakFile + "." + varSuffix + ".vcf"
			cmd_bt_vcf += " -b " + peakFile + ".width" + testSizeArg
			cmd_bt_vcf += " > " + peakFile + ".width" + testSizeArg + "." + varSuffix + ".vcf"
			print (cmd_bt_vcf)
			subprocess.getoutput(cmd_bt_vcf)

			cmd_bt_vcf = bdtArg + " intersect -header -wa -a " + peakFile + ".nonDnD.vcf"
			cmd_bt_vcf += " -b " + peakFile + ".width" + testSizeArg
			cmd_bt_vcf += " > " + peakFile + ".width" + testSizeArg + ".nonDnD.vcf"
			print (cmd_bt_vcf)
			subprocess.getoutput(cmd_bt_vcf)

	if modeArg in ["sea", "homer2"]:
		peak_target_resized = "/motif_TF.narrowPeak.width"
		peak_target_centric = "/motif_TF.narrowPeak.motifCentered.width"
		peak_bkgd_resized = "/motif_bkgd.narrowPeak.width"
		out_resized = "/resizedPeakCounts"
		out_centric = "/motifCenteredCounts"
	elif modeArg == "chip":
		peak_target_resized = "/motif_TF.narrowPeak.width"
		peak_target_centric = "/motif_TF.summits.bed.width"
		peak_bkgd_resized = "/motif_bkgd.summits.bed.width"
		out_resized = "/resizedPeakCounts"
		out_centric = "/summitCenteredCounts"

	if not os.path.exists(outputArg + sampleArg + out_resized):
		os.mkdir(outputArg + sampleArg + out_resized)
	motifCenteredCounts(outputArg + sampleArg + peak_target_resized + testSizeArg, outputArg + sampleArg + peak_target_resized + testSizeArg + "." + varSuffix + ".vcf",
			    outputArg + sampleArg + peak_bkgd_resized + testSizeArg, outputArg + sampleArg + peak_bkgd_resized + testSizeArg  + "." + varSuffix + ".vcf",
			    outputArg + sampleArg + out_resized + "/res", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "Rscript ~/DnD/src/main/tools.R plot_footprint " + outputArg + sampleArg + out_resized
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + sampleArg + out_resized + "_nonDnD"):
		os.mkdir(outputArg + sampleArg + out_resized + "_nonDnD")
	motifCenteredCounts(outputArg + sampleArg + peak_target_resized + testSizeArg, outputArg + sampleArg + peak_target_resized + testSizeArg + ".nonDnD.vcf",
			    outputArg + sampleArg + peak_bkgd_resized + testSizeArg, outputArg + sampleArg + peak_bkgd_resized + testSizeArg  + ".nonDnD.vcf",
			    outputArg + sampleArg + out_resized + "_nonDnD/res", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "Rscript ~/DnD/src/main/tools.R plot_footprint " + outputArg + sampleArg + out_resized + "_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + sampleArg + out_centric):
		os.mkdir(outputArg + sampleArg + out_centric)
	motifCenteredCounts(outputArg + sampleArg + peak_target_centric + testSizeArg, outputArg + sampleArg + peak_target_centric + testSizeArg + "." + varSuffix + ".vcf",
			    outputArg + sampleArg + peak_bkgd_resized + testSizeArg, outputArg + sampleArg + peak_bkgd_resized + testSizeArg  + "." + varSuffix + ".vcf",
			    outputArg + sampleArg + out_centric + "/res", testSizeArg, peakSampleSizeArg, "DnD", varArgs)
	cmd_plt = "Rscript ~/DnD/src/main/tools.R plot_footprint " + outputArg + sampleArg + out_centric
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if not os.path.exists(outputArg + sampleArg + out_centric + "_nonDnD"):
		os.mkdir(outputArg + sampleArg + out_centric + "_nonDnD")
	motifCenteredCounts(outputArg + sampleArg + peak_target_centric + testSizeArg, outputArg + sampleArg + peak_target_centric + testSizeArg + ".nonDnD.vcf",
			    outputArg + sampleArg + peak_bkgd_resized + testSizeArg, outputArg + sampleArg + peak_bkgd_resized + testSizeArg  + ".nonDnD.vcf",
			    outputArg + sampleArg + out_centric + "_nonDnD/res", testSizeArg, peakSampleSizeArg, "nonDnD", varArgs)
	cmd_plt = "Rscript ~/DnD/src/main/tools.R plot_footprint " + outputArg + sampleArg + out_centric + "_nonDnD"
	print (cmd_plt)
	subprocess.getoutput(cmd_plt)

	if modeArg in ["sea", "homer2"]:
		countOnlySnvs(outputArg + sampleArg + "/motif_TF.narrowPeak.motifCentered.width" + testSizeArg + ".vcf", 
				outputArg + sampleArg + "/motif_TF.narrowPeak.motifCentered.width" + testSizeArg, outputArg + sampleArg + "/target")
		countOnlySnvs(outputArg + sampleArg + "/motif_bkgd.narrowPeak.vcf", 
				outputArg + sampleArg + "/motif_bkgd.narrowPeak.width" + testSizeArg, outputArg + sampleArg + "/background")
	elif modeArg == "chip":
		countOnlySnvs(outputArg + sampleArg + "/motif_TF.summits.bed.width" + testSizeArg + ".vcf", 
				outputArg + sampleArg + "/motif_TF.summits.bed.width" + testSizeArg, outputArg + sampleArg + "/target")
		countOnlySnvs(outputArg + sampleArg + "/motif_bkgd.summits.bed.width" + testSizeArg + ".vcf", 
				outputArg + sampleArg + "/motif_bkgd.summits.bed.width" + testSizeArg, outputArg + sampleArg + "/background")
	cmdsnr = "Rscript ~/DnD/src/main/tools.R plot_context"
	cmdsnr += " " + outputArg + sampleArg + "/target.countNormalized.txt"
	cmdsnr += " " + outputArg + sampleArg + "/background.countNormalized.txt"
	cmdsnr += " " + outputArg + sampleArg + "/"
	print (cmdsnr, "\n")
	subprocess.getoutput(cmdsnr)

	### di-/tri-nucleotide context ###
	if modeArg in ["sea", "homer2"]:
		tfVcf = outputArg + sampleArg + "/motif_TF.narrowPeak.motifCentered.width" + testSizeArg + ".vcf"
		bkgdVcf = outputArg + sampleArg + "/motif_bkgd.narrowPeak.width" +testSizeArg + ".vcf" 
	elif modeArg == "chip":
		tfVcf = outputArg + sampleArg + "/motif_TF.summits.bed.width" + testSizeArg + ".vcf"
		bkgdVcf = outputArg + sampleArg + "/motif_bkgd.summits.bed.width" + testSizeArg + ".vcf"

	cmd_v2b_tf = v2bArg + " < " + tfVcf + " > " + tfVcf + ".bed"
	print (cmd_v2b_tf)
	subprocess.getoutput(cmd_v2b_tf)
	cmd_v2b_bkgd = v2bArg + " < " + bkgdVcf + " > " + bkgdVcf + ".bed"
	print (cmd_v2b_bkgd)
	subprocess.getoutput(cmd_v2b_bkgd)

	cmd_ext_bed = "Rscript ~/DnD/src/main/tools.R extend_bed"
	cmd_ext_bed_tf = cmd_ext_bed + " " + tfVcf + ".bed"
	cmd_ext_bed_tf += " " + "5" 
	cmd_ext_bed_tf += " " + varSuffix
	subprocess.getoutput(cmd_ext_bed_tf)
	cmd_ext_bed_bkgd = cmd_ext_bed + " " + bkgdVcf + ".bed"
	cmd_ext_bed_bkgd += " " + "5"
	cmd_ext_bed_bkgd += " " + varSuffix
	subprocess.getoutput(cmd_ext_bed_bkgd)

	target_fas = list()
	background_fas = list()
	for var in varSuffix.split("_"):
		cmd_bdt = bdtArg + " getfasta -fi ~/DnD/genome/genome.fa"
		cmd_bdt_tf = cmd_bdt + " -fo " + tfVcf + "." + var + "_+-5bp.fasta"
		cmd_bdt_tf += " -bed " + tfVcf + "." + var + "_+-5bp.bed"
		print (cmd_bdt_tf)
		subprocess.getoutput(cmd_bdt_tf)
	
		cmd_bdt_bkgd = cmd_bdt + " -fo " + bkgdVcf + "." + var + "_+-5bp.fasta"
		cmd_bdt_bkgd +=  " -bed " + bkgdVcf + "." + var + "_+-5bp.bed"
		print (cmd_bdt_bkgd)
		subprocess.getoutput(cmd_bdt_bkgd)

		if var == "GA": ###
			revComp(tfVcf + "." + var + "_+-5bp.fasta")
			revComp(bkgdVcf + "." + var + "_+-5bp.fasta")
			target_fas.append(tfVcf + "." + var + "_+-5bp" + "_revComp.fasta")
			background_fas.append(bkgdVcf + "." + var + "_+-5bp" + "_revComp.fasta")
		else:
			target_fas.append(tfVcf + "." + var + "_+-5bp.fasta")
			background_fas.append(bkgdVcf + "." + var + "_+-5bp.fasta")

	if not os.path.exists(outputArg + sampleArg + "/context_di"):
		os.mkdir(outputArg + sampleArg + "/context_di")
	getNt(target_fas, outputArg + sampleArg + "/context_di/target_di.txt", 5, 1, 0)
	getNt(background_fas, outputArg + sampleArg + "/context_di/background_di.txt", 5, 1, 0)

	if not os.path.exists(outputArg + sampleArg + "/context_tri"):
		os.mkdir(outputArg + sampleArg + "/context_tri")
	getNt(target_fas, outputArg + sampleArg + "/context_tri/target_tri.txt", 5, 2, 0)
	getNt(background_fas, outputArg + sampleArg + "/context_tri/background_tri.txt", 5, 2, 0)

	cmd_plt = "Rscript ~/DnD/src/main/tools.R ntcontext"
	cmd_plt_tf_di = cmd_plt + " " + outputArg + sampleArg + "/context_di/target_di.txt 2"
	subprocess.getoutput(cmd_plt_tf_di)
	cmd_plt_bkgd_di = cmd_plt + " " + outputArg + sampleArg + "/context_di/background_di.txt 2"
	subprocess.getoutput(cmd_plt_bkgd_di)
	cmd_plt_tf_tri = cmd_plt + " " + outputArg + sampleArg + "/context_tri/target_tri.txt 3"
	subprocess.getoutput(cmd_plt_tf_tri)
	cmd_plt_bkgd_tri = cmd_plt + " " + outputArg + sampleArg + "/context_tri/background_tri.txt 3"
	subprocess.getoutput(cmd_plt_bkgd_tri)



