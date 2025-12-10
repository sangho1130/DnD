
import argparse
import os
import sys
import subprocess
import gzip
import random
from multiprocessing import Process
from main.dnd_acc import *


def picard_markDup(pcdArg, dirArg, outputArg, bamArg):
	cmd = pcdArg + " MarkDuplicates"
	cmd += ' I=' + dirArg + bamArg
	cmd += ' O=' + outputArg + "__tmp_rmDup__" + bamArg.lower()
	cmd += ' M=' + outputArg + "step1_picard_deduplication.txt" 
	cmd += ' REMOVE_DUPLICATES=true'
	return cmd, "__tmp_rmDup__" + bamArg.lower()


def samtools_flt(smtArg, outputArg, bamArg, threadArg, mapqArg, chromArg, otherArg, seArg):
	cmd = smtArg + " view -b"
	cmd += ' --threads ' + threadArg
	cmd += ' -q ' + mapqArg # 20 by default
	if not seArg:
		cmd += ' -f 0x2' # read mapped in proper pair
	cmd += ' -F 256' # leave only primary alignment
	if not chromArg:
		cmd += ' -L ' + '~/DnD/genome/GRCh38.p14.genome.chrs.bed'
	if otherArg:
		cmd += ' ' + otherArg
	cmd += ' ' + outputArg + bamArg
	cmd += ' > ' + outputArg + bamArg.replace('__tmp_rmDup__', '')
	return cmd, bamArg.replace('__tmp_rmDup__', '')


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



