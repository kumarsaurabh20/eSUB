#!/usr/bin/env python2
from __future__ import (print_function, division)
from datetime import datetime
from Bio import SeqIO
from FastaIndex import FastaIndex
import sys, os, gzip, math, subprocess

__author__ = Kumar
__desc__ = "Align genomes or dna segments using HMMER/BLAT"

def execute_hmmer():
	pass

def execute_blat():
	pass

def mask_dna():
	pass

def merge_fasta():
	pass

def main():
	import argparse	

	usage   = "%(prog)s -v" #usage=usage, 
	parser  = argparse.ArgumentParser(description=desc, epilog=epilog, formatter_class=argparse.RawTextHelpFormatter)
	####### 	 
    	parser.add_argument('--version', action='version', version='1.01d')   
	parser.add_argument("-v", "--verbose", default=False, action="store_true", help="verbose")    
	parser.add_argument("-i", "-f", "--fasta", nargs="+", type=file, help="FASTA file(s)")
	parser.add_argument("-t", "--threads", default=4, type=int, help="max threads to run [%(default)s]")
	parser.add_argument("--identity",    default=0.51, type=float, help="min. identity   [%(default)s]")
	parser.add_argument("--overlap",     default=0.66, type=float, help="min. overlap    [%(default)s]")
	parser.add_argument("--joinOverlap", default=200, type=int, help="min. end overlap to join two contigs [%(default)s]")
	parser.add_argument("--endTrimming", default=33, type=int, help="max. end trim on contig join [%(default)s]")
	parser.add_argument("--minLength",   default=200, type=int, help="min. contig length [%(default)s]")
	#######    
	o = parser.parse_args()
	if o.verbose:
        	sys.stderr.write("Options: %s\n"%str(o))	

if __name__ =='__main__':
	t0 = datetime.now()
	try:	
        	main()
	except KeyboardInterrupt:
        	sys.stderr.write("\nCtrl-C pressed!      \n")

	dt = datetime.now()-t0
	sys.stderr.write("#Time elapsed: %s\n"%dt)
