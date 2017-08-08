#!/usr/bin/env python2
from __future__ import print_function
import pysam
import pandas as pd
from collections import OrderedDict
from operator import itemgetter

bamfile = pysam.AlignmentFile("pt_with_control_data.sorted.bam", "rb")

if bamfile.check_index():
    print("Samfile is indexed!!")
else:
    print("index unavailable. Create index for your bam file!!")

fasta_map = {}
with open("fasta_length.tab", "r") as fl:
    for line in fl:
        temp = line.split("\t")
        fasta_map[temp[0].strip()] = temp[1].strip()

count = OrderedDict()
window_size = 19644

binCov = 0
#temp_count = {}
#m = 0
#n = 0
# for each in fasta_map.keys():
# 	temp_count = {}
# 	temp_key = ""
# 	startPos = 1
# 	stopPos = window_size
# 	#print("%s ==> %s"%(each, fasta_map[each]))
# 	length = int(fasta_map[each])
# 	if int(length) <= int(window_size):
# 		#print("%s length is shorter than window length"%(each))
# 		#m = m + 1
# 		binCov = 0
# 		for cov in bamfile.fetch(each, startPos, length):
# 			binCov = binCov + 1
# 			#print("%s\t%d\t%d\t%d"%(each,startPos,final_pos,binCov))
# 		temp_key = "%s.%s%d" % (each, "bin", 1)
# 		temp_count[temp_key] = binCov
# 		count[each] = temp_count
# 		#print("################## Printing filled temp")
# 		#print(temp_count)
# 		#temp_count.clear()
# 		#print("&&&&&&&&&&&&&&&&&& Printing empty temp")
# 		#print(temp_count)
#
# 	elif int(length) > int(window_size):
# 		#print("Reading scaffold: %s"%(each))
# 		#n = n + 1
# 		binCov = 0
# 		numOfBins = int(length / window_size)
# 		#print(numOfBins)
# 		penumPos = numOfBins * window_size
# 		for i in xrange(1, numOfBins + 1):
# 			#penumPos = numOfBins * window_size
# 			if i <= numOfBins:
# 				binCov = 0
# 				for cov in bamfile.fetch(each, startPos, stopPos):
# 					binCov = binCov + 1
# 				temp_key = "%s.%s%d"%(each,"bin", i)
# 				temp_count[temp_key] = binCov
# 				#print("Inserting binCov %d"%binCov)
# 				#print(temp_count)
# 			elif i > numOfBins:
# 				binCov = 0
# 				startPos = penumPos + 1
# 				for cov in bamfile.fetch(each, startPos, length):
# 					binCov = binCov + 1
# 				temp_key = "%s.%s%d" % (each, "bin", i)
# 				temp_count[temp_key] = binCov
# 			else:
# 				continue
# 			startPos = window_size + 1
# 			stopPos = startPos + window_size
# 		count[each] = temp_count
# 		#temp_count.clear()
# 	else:
# 		continue
# 	#temp_count.clear()

step = int(window_size / 2)
for each in fasta_map.keys():
    temp_count = OrderedDict()
    temp_key = ""
    startPos = 1
    stopPos = window_size
    length = int(fasta_map[each])
    if int(length) <= int(window_size):
        binCov = 0
        for cov in bamfile.fetch(each, startPos, length):
            binCov = binCov + 1
        temp_key = "%s.%s%d" % (each, "bin", 1)
        data = [startPos, length, binCov]
        temp_count[temp_key] = data
        count[each] = temp_count

    elif int(length) > int(window_size):
        # print("Reading scaffold: %s"%(each))
        # n = n + 1
        binCov = 0
        numOfBins = int(length / step)
        # print(numOfBins)
        penumPos = numOfBins * step
        for i in xrange(1, (numOfBins)+2):
            # penumPos = numOfBins * window_size
            if i == 1:
                binCov = 0
                for cov in bamfile.fetch(each, startPos, stopPos):
                    binCov = binCov + 1
                temp_key = "%s.%s%d" % (each, "bin", i)
                data = [startPos, stopPos, binCov]
                temp_count[temp_key] = data
                #print(temp_count)
            elif i <= numOfBins:
                binCov = 0
                for cov in bamfile.fetch(each, startPos, stopPos):
                    binCov = binCov + 1
                temp_key = "%s.%s%d" % (each, "bin", i)
                data = [startPos, stopPos, binCov]
                temp_count[temp_key] = data
                #print(temp_count)
            elif i > numOfBins:
                binCov = 0
                startPos = penumPos + 1
                stopPos = length
                for cov in bamfile.fetch(each, startPos, length):
                    binCov = binCov + 1
                temp_key = "%s.%s%d" % (each, "bin", i)
                data = [startPos, stopPos, binCov]
                temp_count[temp_key] = data
                #print(temp_count)
            else:
                continue

            startPos = (stopPos - step) + 1
            stopPos = (startPos + window_size) - 1
        count[each] = temp_count
    # temp_count.clear()
    else:
        continue
        # temp_count.clear()

ch = open("data.cov", "w")
for key, value in count.items():
    for k,v in value.items():
        print("%s\t%d\t%d\t%d"%(key, v[0], v[1], v[2]))
        if ch.tell() == 0:
            ch.write("%s\t%s\t%s\t%s\n"%("Chromosome", "Start", "End", "Coverage"))
            ch.write("%s\t%d\t%d\t%d\n" % (key, v[0], v[1], v[2]))
            #continue
        else:
            ch.write("%s\t%d\t%d\t%d\n"%(key, v[0], v[1], v[2]))
ch.close()

bamfile.close()
