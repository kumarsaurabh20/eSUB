#!/usr/bin/env python2
from __future__ import print_function
import pysam
from scipy.stats import norm
from collections import OrderedDict
import sys, os
import subprocess


##methods
def get_T_for_max_ws(p_value):
    return norm.ppf(1 - 0.5 * p_value)

def get_T_for_min_ws(p_value):
    return norm.ppf(0.5 * p_value)

def get_r_for_max_ws(r):
    return r**2

def get_r_for_min_ws(r):
    return (1 / r**2)

def get_max_window_size(Nx, Ny, Tmax,rmax, G):
    return (( Nx * rmax**2 + Ny ) * G * Tmax**2 ) / ((1 - rmax)**2 * Nx * Ny)

def get_min_window_size(Nx, Ny, Tmin, rmin, G):
    return (( Nx * rmin**2 + Ny ) * G * Tmin**2 ) / ((1 - rmin)**2 * Nx * Ny)

def check_alignment_file(bam):
    bamfile = pysam.AlignmentFile(bam, "rb")

    if bamfile.check_index():
        print("bamfile is indexed!")
    else:
        print("index unavailable for bamfile. Create index and try again!!")

    return bamfile

def create_fasta_map(ref_fasta):
    fasta_map = {}

    with open("fasta_length.tab", "r") as fl:
        for line in fl:
            temp = line.split("\t")
            fasta_map[temp[0].strip()] = temp[1].strip()
    return fasta_map

def count_non_sliding(fasta_map, bamHandler, window_size):
    count_map = {}
    for each in fasta_map.keys():
        temp_count = {}
        temp_key = ""
        startPos = 1
        stopPos = window_size
        # print("%s ==> %s"%(each, fasta_map[each]))
        length = int(fasta_map[each])
        if int(length) <= int(window_size):
            # print("%s length is shorter than window length"%(each))
            # m = m + 1
            binCov = 0
            for cov in bamHandler.fetch(each, startPos, length):
                binCov = binCov + 1
            # print("%s\t%d\t%d\t%d"%(each,startPos,final_pos,binCov))
            temp_key = "%s.%s%d" % (each, "bin", 1)
            temp_count[temp_key] = binCov
            count_map[each] = temp_count

        elif int(length) > int(window_size):
            print("Reading scaffold: %s" % (each))
            # n = n + 1
            binCov = 0
            numOfBins = int(length / window_size)
            print(numOfBins)
            penumPos = numOfBins * window_size

            for i in xrange(1, numOfBins + 1):
                # penumPos = numOfBins * window_size
                if i <= numOfBins:
                    binCov = 0
                    for cov in bamHandler.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    temp_count[temp_key] = binCov
                    print("Inserting binCov %d" % binCov)
                    print(temp_count)
                elif i > numOfBins:
                    binCov = 0
                    startPos = penumPos + 1
                    for cov in bamHandler.fetch(each, startPos, length):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    temp_count[temp_key] = binCov
                else:
                    continue
                startPos = window_size + 1
                stopPos = startPos + window_size
            count_map[each] = temp_count

        else:
            continue
    return count_map

def count_sliding(fasta_map, bamHandler, window_size, step):
    count_map = OrderedDict()
    for each in fasta_map.keys():
        temp_count = OrderedDict()
        temp_key = ""
        startPos = 1
        stopPos = window_size
        # print("%s ==> %s"%(each, fasta_map[each]))
        length = int(fasta_map[each])
        if int(length) <= int(window_size):
            # print("%s length is shorter than window length"%(each))
            binCov = 0
            for cov in bamHandler.fetch(each, startPos, length):
                binCov = binCov + 1
            # print("%s\t%d\t%d\t%d"%(each,startPos,final_pos,binCov))
            temp_key = "%s.%s%d" % (each, "bin", 1)
            data = [startPos, length, binCov]
            temp_count[temp_key] = data
            # print(temp_count)
            count_map[each] = temp_count

        elif int(length) > int(window_size):
            # print("Reading scaffold: %s"%(each))
            # n = n + 1
            binCov = 0
            numOfBins = int(length / step)
            # print(numOfBins)
            penumPos = numOfBins * step
            for i in xrange(1, (numOfBins) + 2):
                # penumPos = numOfBins * window_size
                if i == 1:
                    binCov = 0
                    for cov in bamHandler.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                elif i <= numOfBins:
                    binCov = 0
                    for cov in bamHandler.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                elif i > numOfBins:
                    binCov = 0
                    startPos = penumPos + 1
                    stopPos = length
                    for cov in bamHandler.fetch(each, startPos, length):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                else:
                    continue

                startPos = (stopPos - step) + 1
                stopPos = (startPos + window_size) - 1
            count_map[each] = temp_count
        # temp_count.clear()
        else:
            continue
    return count_map       # temp_count.clear()

def get_windows(refBam, qryBam, p_value, log2, alpha, genome_size):
    # total_test is the total number of reads in test set
    # total_ref is the total number of reads in test_ref
    total_query_reads = refBam.count(until_eof=True)
    total_ref_reads = qryBam.count(until_eof=True)
    #
    r_bt = get_T_for_max_ws(p_value)
    r_st = get_T_for_min_ws(p_value)
    log2 = abs(log2)
    brp = get_r_for_max_ws(log2)
    srp = get_r_for_min_ws(log2)
    bw = get_max_window_size(total_query_reads, total_ref_reads, r_bt, brp, genome_size)
    sw = get_min_window_size(total_query_reads, total_ref_reads, r_st, srp, genome_size)
    #
    window_size = max(bw, sw)
    window_size *= alpha
    window_size = '{:.0}'.format(window_size)
    # window_size = 19644
    step = int(window_size / 2)
    #
    return window_size, step

def write_cov(count_map={}, filename=""):
    ch = open(filename, "w")
    for key, value in count_map.items():
        for k, v in value.items():
            print("%s\t%d\t%d\t%d" % (key, v[0], v[1], v[2]))
            if ch.tell() == 0:
                ch.write("%s\t%s\t%s\t%s\n" % ("Chromosome", "Start", "End", "Coverage"))
                ch.write("%s\t%d\t%d\t%d\n" % (key, v[0], v[1], v[2]))
            else:
                ch.write("%s\t%d\t%d\t%d\n" % (key, v[0], v[1], v[2]))
    ch.close()
    print("Coverage file is written successfully!!")

def normalize_counts(refBam, qryBam):
    refCounts = refBam.count(until_eof=True)
    qryCounts = qryBam.count(until_eof=True)
    factor = refCounts / qryCounts
    return factor

#https://stackoverflow.com/questions/16932825/why-non-default-arguments-cant-follows-default-argument
def count_coverage(refBam, qryBam, ref_fasta, p_value, log2, alpha, genome_size, filename, steps= False):
    #
    subtraction_map = {}
    refBamHandler = check_alignment_file(refBam)
    qryBamHandler = check_alignment_file(qryBam)
    #
    fasta_map = create_fasta_map(ref_fasta)
    window_size, steps= get_windows(p_value, log2, alpha, genome_size)
    # print(samfile.count(until_eof=True))
    binCov = 0
    if steps == True:
        count_ref_map = count_sliding(fasta_map, refBamHandler, window_size, steps)
        count_qry_map = count_sliding(fasta_map, qryBamHandler, window_size, steps)
    else:
        count_ref_map = count_non_sliding(fasta_map, refBamHandler, window_size)
        count_qry_map = count_non_sliding(fasta_map, qryBamHandler, window_size)
    #
    write_cov(subtraction_map, filename)
    refBamHandler.close()
    qryBamHandler.close()


#https://github.com/JohnLonginotto/pybam

#pysam.sort("-o", "output.bam", "ex1.bam")
#samtools sort -o output.bam ex1.bam
#pysam.sort("-m", "1000000", "-o", "output.bam", "ex1.bam")
# tabixfile = pysam.TabixFile("example.gtf.gz")
# for gtf in tabixfile.fetch("chr1", 1000, 2000):
#     print (gtf.contig, gtf.start, gtf.end, gtf.gene_id)