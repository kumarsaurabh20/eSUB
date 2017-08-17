#!/usr/bin/env python2
from __future__ import print_function
import sys
import subprocess
import os
import argparse
import shutil
import datetime
import ConfigParser

def main():
	startTime = datetime.datetime.now()
	global args
	args = vars(getArguments())
	checkArguments()
	printStartMessage()
	readConfigFile()
	makeOutputDirectory()

def makeOutputDirectory():
	fullOutPath = os.getcwd() + '/' + args['o'] + '/'

	global outDir
	outDir = os.path.dirname(fullOutPath)

	if not os.path.exists(outDir):
		os.makedirs(outDir)

def deleteTemporaryDirectories(iterDir):
	indexDir = iterDir + '/1_mapping_index'
	pairedDir = iterDir + '/2-paired_read_alignments'
	unpairedDir = iterDir + '/2-unpaired_read_alignments'
	assemblyDir = iterDir + '/3-assembly'

	shutil.rmtree(indexDir)
	if os.path.exists(pairedDir):
		shutil.rmtree(pairedDir)
	if os.path.exists(unpairedDir):
		shutil.rmtree(unpairedDir)
		shutil.rmtree(assemblyDir)

# the program.
if __name__ == '__main__':
	main()















