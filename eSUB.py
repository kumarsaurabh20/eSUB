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

def getArguments():

	parser = argparse.ArgumentParser(description=printStartMessage(), add_help=False)

	required = parser.add_argument_group('Required arguments')
	optional = parser.add_argument_group('Optional arguments')

	reads = parser.add_argument_group('Read arguments', 'Read files must be given as interleave PE files. r1 = Reference and r2 = Query')
	contigs = parser.add_argument_group('Assembled sequence arguments', 'Contigs as multiple fasta files. c1 = Reference and c2 = Query')
	aligns = parser.add_argument_group('Alignment file argumetns', 'Alignment file in sorted BAM format. a1 = Reference and a2 = Query')

	optional.add_argument("-h", "--help", action="help", help="show this help message and exit")
	optional.add_argument("-k", "--keep", help="keep all temporary files", action="store_true")
	optional.add_argument("-f", "--fraction", help="Kmer fraction between 0 and 1. Lets you decide how much kmer fraction is required to decide read uniqueness")
	optional.add_argument("-c", "--config", help="Configuration file which specifies commands")

	required.add_argument("-m", metavar="mode", help="Mode specifies which method should be used with the respective datasets. 1 = k-mer based subtraction of read dataset 2 = k-mer based subtraction of sequence dataset 3 = Similarity based subtraction of sequence dataset 4 = Alignment file subtraction", required=True)
	required.add_argument("-o", metavar="OUTDIR", help="the output directory", required=True)
	required.add_argument("-p", metavar="PREFIX", help="all the output file is annoted with the user specified prefix", required=True)

	reads.add_argument("-r1", metavar="FIRST", help="Reference interleaved reads file")
	reads.add_argument("-r2", metavar="SECOND", help="Query interleaved reference file")

	contigs.add_argument("-s1", metavar="FIRST", help="Reference sequences file")
	contigs.add_argument("-s2", metavar="SECOND", help="Query sequences file")

	aligns.add_argument("-a1", metavar="FIRST", help="Reference alignment file")
	aligns.add_argument("-a2", metavar="SECOND", help="Query alignment file")

	return parser.parse_args()

def checkArguments():
	global args

	if (args['r1'] == None and args['r2'] != None):
		print('If a query dataset is given, then a reference dataset is also required.')
		exit()
	if (args['r1'] != None and args['r2'] == None):
		print('If a first mate file is given, then a second mate file is also required.')
		exit()
	if args['r1'] != None and not os.path.isfile(args['r1']):
		print('The reference dataset could not be found.')
		exit()
	if args['r2'] != None and not os.path.isfile(args['r2']):
		print('The query dataset could not be found.')
		exit()

	if (args['c1'] == None and args['c2'] != None):
		print('If a query dataset is given, then a reference dataset is also required.')
		exit()
	if (args['c1'] != None and args['c2'] == None):
		print('If a first mate file is given, then a second mate file is also required.')
		exit()
	if args['c1'] != None and not os.path.isfile(args['c1']):
		print('The reference dataset could not be found.')
		exit()
	if args['c2'] != None and not os.path.isfile(args['c2']):
		print('The query dataset could not be found.')
		exit()

	if (args['a1'] == None and args['a2'] != None):
		print('If a query dataset is given, then a reference dataset is also required.')
		exit()
	if (args['a1'] != None and args['a2'] == None):
		print('If a first mate file is given, then a second mate file is also required.')
		exit()
	if args['a1'] != None and not os.path.isfile(args['a1']):
		print('The reference dataset could not be found.')
		exit()
	if args['a2'] != None and not os.path.isfile(args['a2']):
		print('The query dataset could not be found.')
		exit()


def readConfigFile():

	global args
	configFileName = args['c']
	if not os.path.isfile(configFileName):
		print("\nERROR: The configuration file could not be found.\n")
		exit()

	config = ConfigParser.ConfigParser()
	config.read(configFileName)

	if config.has_option('Reads', 'reference'):
        	args['r1'] = config.get('Mapping', 'reference').strip()
	if config.has_option('Reads', 'query'):
		args['r2'] = config.get('Mapping', 'query').strip()
	if config.has_option('Reads','mode'):
		args['m'] = config.get('Reads','mode').strip()

	if config.has_option('Contigs', 'reference'):
		args['c1'] = config.get('Contigs', 'reference').strip()
	if config.has_option('Contigs', 'query'):
		args['c2'] = config.get('Contigs', 'query').strip()
	if config.has_option('Contigs','mode'):
		args['m'] = config.get('Contigs','mode').strip()

	if config.has_option('Alignment', 'reference'):
		args['a1'] = config.get('Alignment', 'reference').strip()
	if config.has_option('Alignment', 'query'):
		args['a2'] = config.get('Alignment', 'query').strip()
	if config.has_option('Alignment','mode'):
		args['m'] = config.get('Alignment','mode').strip()

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

def printStartMessage():
	print(
	"""
	######################################################                            
		           _____  _    _  ____  
		          / ____|| |  | ||  _ \ 
		     ___ | (___  | |  | || |_) |
		    / _ \ \___ \ | |  | ||  _ < 
		   |  __/ ____) || |__| || |_) |
		    \___||_____/  \____/ |____/ 

	######################################################                                	                                
	"""
	)
	print("")
	print("eSUB::Electronic Subtraction of Read Data Using K-mers")
	print("Kumar Saurabh Singh (k.saurabh-singh@exeter.ac.uk)")
	print("")

def convertTimeDeltaToReadableString(timeDelta):
	seconds = timeDelta.seconds
	hours = timeDelta.days * 24
	hours += seconds // 3600
	seconds = seconds % 3600
	minutes = seconds // 60
	seconds = seconds % 60
	seconds += timeDelta.microseconds / 1000000.0
	secondString = "{:.1f}".format(seconds)

	returnString = ""
	if hours > 0:
		return str(hours) + ' h, ' + str(minutes) + ' min, ' + secondString + ' s'
	if minutes > 0:
		return str(minutes) + ' min, ' + secondString + ' s'
	return secondString + ' s'

def printEndMessage():
	print("############################################")
	print("eSUB::Read data subtraction is finished now!")
	print("Total time to complete: %d"% convertTimeDeltaToReadableString())
	print("############################################")

def getDateTimeString():
	return datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

# Standard boilerplate to call the main() function to begin
# the program.
if __name__ == '__main__':
	main()















