"""
This module subtract two bam files by counting coverage using sliding/non-sliding windows on the reference and query genomes
and returning P-values for regions which are significantly different in the query genome relative to reference.
"""
from __future__ import print_function
from Bio import SeqIO
from collections import OrderedDict
import pysam, sys, os
import numpy as np
import pandas as pd
import scipy.stats

# pd.set_option('display.height', 1000)
# pd.set_option('display.max_rows', 500)
# pd.set_option('display.max_columns', 500)
# pd.set_option('display.width', 1000)

__all__ = ('subtractBAMs',)

#bam sanity checks
def check_alignment_file(bam):

    if not os.path.isfile(bam):
        print("\nERROR: The bam file could not be found.\n")
        exit()

    try:
        bamfile = pysam.AlignmentFile(bam, "rb")

        if bamfile.check_index():
            print("bamfile is indexed!")
        else:
            print("index unavailable for bamfile. Create index and try again!!")
            exit()
    except IOError:
        print("Program encountered error while with bam file!!")
    return bamfile

#fasta sanity checks
def is_fasta(filename):

    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

#Create fasta dictionary
def create_fasta_map(ref_fasta):

    assert is_fasta(ref_fasta), "Check the format of reference fasta file!!"
    fasta_map = {}
    try:
        for s in SeqIO.parse(ref_fasta, "fasta"):
            # temp = s.id.split("")[0]
            temp = str(s.id)
            fasta_map[temp] = len(s.seq)
            #print("%s\t%i" % (temp, len(s.seq)))
    except IOError:
        print("Program encountered error while creating fasta dictionary. Check the reference fasta file!")

    return fasta_map


##calculate algo values
def get_T_for_max_ws(p_value):
    #assert isinstance(p_value, float), "input value is not floating point!"
    return scipy.stats.norm.ppf(1 - 0.5 * p_value)

def get_T_for_min_ws(p_value):
    #assert isinstance(p_value, float), "input value is not floating point!"
    return scipy.stats.norm.ppf(0.5 * p_value)

def get_r_for_max_ws(r):
    #assert r, "Check the log2 input value!!"
    return 2**r

def get_r_for_min_ws(r):
    #assert r, "Check the log2 input value!!"
    return (1 / 2**r)

def get_max_window_size(Nx, Ny, Tmax,rmax, G):
    try:
        bw = (( Nx * rmax**2 + Ny ) * G * Tmax**2 ) / ((1 - rmax)**2 * Nx * Ny)
    except:
        print("Program encountered error while calculating max window size!!")
    return bw

def get_min_window_size(Nx, Ny, Tmin, rmin, G):
    try:
        sw = (( Nx * rmin**2 + Ny ) * G * Tmin**2 ) / ((1 - rmin)**2 * Nx * Ny)
    except:
        print("Program encountered error while calculating min window size!!")
    return sw

#estimate the window size using input log2 and p-value
def get_windows(refBam, qryBam, p_value, log2, alpha, genome_size):
    # total_test is the total number of reads in test set
    # total_ref is the total number of reads in test_ref
    total_query_reads = refBam.count(until_eof=True)
    #print(total_query_reads)
    total_ref_reads = qryBam.count(until_eof=True)
    #print(total_ref_reads)

    r_bt = get_T_for_max_ws(p_value)
    #print(r_bt)
    r_st = get_T_for_min_ws(p_value)
    log2 = abs(log2)
    brp = get_r_for_max_ws(log2)
    #print(brp)
    srp = get_r_for_min_ws(log2)
    bw = get_max_window_size(total_query_reads, total_ref_reads, r_bt, brp, genome_size)
    sw = get_min_window_size(total_query_reads, total_ref_reads, r_st, srp, genome_size)
    #print(bw)
    #
    window_size = max(bw, sw)
    window_size *= alpha
    window_size = float(window_size)
    # window_size = 19644
    #step = int(window_size / 2)
    #
    return window_size

#calculate coverage using a bed file
def bed_coverage(bamfile, bedfile, prefix="out"):
    keys, start, end, coverage = ([] for i in range(4))
    bamHandler = check_alignment_file(bamfile)
    counter = 0
    outfile="%s_cov.bed"%(prefix)
    out = open(outfile, "a")
    num_lines = sum(1 for line in open(bedfile))
    with open(bedfile, "r") as f:
        for line in f:
            temp = line.strip().split()
            #print(temp)
            startPos = int(temp[1])
            stopPos = int(temp[2])
            scaffold = str(temp[0])
            binCov = 0
            for cov in bamHandler.fetch(scaffold, startPos, stopPos):
                binCov = binCov + 1
            header="chromosome\tStart\tEnd\tCoverage\n"
            if out.tell() == 0:
                out.write(header)
                out.write("%s\t%d\t%d\t%d\n"%(scaffold,startPos,stopPos,binCov))
            else:
                out.write("%s\t%d\t%d\t%d\n"%(scaffold,startPos,stopPos,binCov))
            counter += 1
            #print(os.path.getsize(file_name)/1024+'KB / '+size+' KB downloaded!', end='\r')
            print("Finished %i interval"%(counter), end='\r')
            sys.stdout.flush()
    print("Done!")
    out.close

#calculate coverage using non sliding/overlapping windows
def non_sliding_coverage(bamfile, fasta_map, window_size):
    count_ns = OrderedDict()
    binCov = 0
    for each in fasta_map.keys():
        temp_count = OrderedDict()
        temp_key = ""
        startPos = 1
        stopPos = int(window_size)
        # print("%s ==> %s"%(each, fasta_map[each]))
        length = int(fasta_map[each])
        if int(length) <= int(window_size):
            # print("%s length is shorter than window length"%(each))
            # m = m + 1
            binCov = 0
            for cov in bamfile.fetch(each, startPos, length):
                binCov = binCov + 1
            # print("%s\t%d\t%d\t%d"%(each,startPos,final_pos,binCov))
            temp_key = "%s.%s%d" % (each, "bin", 1)
            data = [startPos, length, binCov]
            temp_count[temp_key] = data
            count_ns[each] = temp_count

        elif int(length) > int(window_size):
            # print("Reading scaffold: %s"%(each))
            # n = n + 1
            binCov = 0
            numOfBins = int(length / window_size)
            # print(numOfBins)
            penumPos = numOfBins * window_size
            for i in xrange(1, numOfBins + 2):
                # penumPos = numOfBins * window_size
                if i <= numOfBins:
                    binCov = 0
                    for cov in bamfile.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, length, binCov]
                    temp_count[temp_key] = data
                # print("Inserting binCov %d"%binCov)
                # print(temp_count)
                elif i > numOfBins:
                    binCov = 0
                    startPos = penumPos + 1
                    for cov in bamfile.fetch(each, startPos, length):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, length, binCov]
                    temp_count[temp_key] = data
                else:
                    continue
                startPos = window_size + 1
                stopPos = startPos + window_size
            count_ns[each] = temp_count
        # temp_count.clear()
        else:
            continue
            # temp_count.clear()
    return count_ns

#calculate coverage using sliding/overlapping windows
def sliding_coverage(bamfile, fasta_map, window_size):
    count_ns = OrderedDict()
    #bamfile = pysam.AlignmentFile(file, "rb")
    #total_reads = bamfile.count(until_eof=True)
    #binCov = 0
    #print(type(window_size))
    #print(window_size)
    print(int(window_size))
    step = int(window_size / 2)
    for each in fasta_map.keys():
        temp_count = OrderedDict()
        temp_key = ""
        startPos = 1
        stopPos = int(window_size)
        # print("%s ==> %s"%(each, fasta_map[each]))
        length = int(fasta_map[each])
        if int(length) <= int(window_size):
            # print("%s length is shorter than window length"%(each))
            binCov = 0
            for cov in bamfile.fetch(each, startPos, length):
                binCov = binCov + 1
            # print("%s\t%d\t%d\t%d"%(each,startPos,final_pos,binCov))
            temp_key = "%s.%s%d" % (each, "bin", 1)
            data = [startPos, length, binCov]
            temp_count[temp_key] = data
            # print(temp_count)
            count_ns[each] = temp_count

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
                    for cov in bamfile.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                elif i <= numOfBins:
                    binCov = 0
                    for cov in bamfile.fetch(each, startPos, stopPos):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                elif i > numOfBins:
                    binCov = 0
                    startPos = penumPos + 1
                    stopPos = length
                    for cov in bamfile.fetch(each, startPos, length):
                        binCov = binCov + 1
                    temp_key = "%s.%s%d" % (each, "bin", i)
                    data = [startPos, stopPos, binCov]
                    temp_count[temp_key] = data
                    # print(temp_count)
                else:
                    continue

                startPos = (stopPos - step) + 1
                stopPos = (startPos + int(window_size) - 1)
            count_ns[each] = temp_count
        # temp_count.clear()
        else:
            continue
            # temp_count.clear()
    return count_ns
    #bamfile.close()

#convert dictionary data in to numpy arrays
def extractData(hash_map):
    keys, start, end, coverage = ([] for i in range(4))
    #ch = open(filename,"a")
    for key, value in hash_map.items():
        for k,v in value.items():
            keys.append(key)
            start.append(v[0])
            end.append(v[1])
            coverage.append(v[2])
    return np.array(keys), np.array(start), np.array(end), np.array(coverage)

#modular function for coverage count
#https://stackoverflow.com/questions/16932825/why-non-default-arguments-cant-follows-default-argument
def count_coverage(refBam, qryBam, ref_fasta, p_value, log2, alpha, genome_size, steps=False):
    #
    ref_keys, qry_keys, ref_start, qry_start, ref_end, qry_end, ref_cov, qry_cov = ([] for i in range(8))
    #
    refBamHandler = check_alignment_file(refBam)
    qryBamHandler = check_alignment_file(qryBam)
    #
    fasta_map = create_fasta_map(ref_fasta)
    window_size = get_windows(refBamHandler, qryBamHandler, p_value, log2, alpha, genome_size)
    # print(samfile.count(until_eof=True))
    binCov = 0
    if steps == True:
        count_ref_map = sliding_coverage(refBamHandler, fasta_map, window_size)
        count_qry_map = sliding_coverage(qryBamHandler, fasta_map, window_size)
        ref_keys, ref_start, ref_end, ref_cov = extractData(count_ref_map)
        qry_keys, qry_start, qry_end, qry_cov = extractData(count_qry_map)
    else:
        count_ref_map = non_sliding_coverage(refBamHandler, fasta_map, window_size)
        count_qry_map = non_sliding_coverage(qryBamHandler, fasta_map, window_size)
        ref_keys, ref_start, ref_end, ref_cov = extractData(count_ref_map)
        qry_keys, qry_start, qry_end, qry_cov = extractData(count_qry_map)
    #
    refBamHandler.close()
    qryBamHandler.close()
    return ref_keys, ref_start, ref_end, qry_cov, ref_cov

#create a panda dataframe
def createCountData(ref_keys, ref_start, ref_end, qry_cov, ref_cov):
    cols = ['chromosome', 'start', 'end', 'test', 'ref']
    data = pd.DataFrame({'chromosome': ref_keys, 'start': ref_start, 'end': ref_end, 'test': qry_cov, 'ref': ref_cov},
                        columns=cols)
    return data

#estimate log2 values and p-values for individual scaffold fragments
def estimateSignificance(data, chr_normalization=False):
    #data = pd.read_table(count_file, header=0)
    print(data.shape)
    data['Position'] = np.floor((data['end'] + data['start']) / 2)
    data['log2'] = np.nan
    data['P_value'] = np.nan
    #print(data)
    chrom = data.chromosome.unique()
    #chr_normalization = False

    if chr_normalization == True:
        for c in chrom:
            sub = data[data['chromosome'] == c]
            test_mean = sub['test'].mean()
            ref_mean = sub['ref'].mean()
            local_norm = test_mean / ref_mean
            local_ratio = sub['test'] / sub['ref']
            data.iloc[sub.index, 6] = np.log2(sub['test'] / sub['ref'] / local_norm)
            positive_sub = data[data['log2'] >= 0]
            #positive_sub = sub.loc[sub['log2'] >= 0, :]
            #print(positive_sub.shape)
            negative_sub = data[data['log2'] < 0]
            #negative_sub = sub.loc[sub['log2'] < 0, :]
            z = positive_sub['test'] / positive_sub['ref']
            zn = negative_sub['test'] / negative_sub['ref']
            # p = 1 - scipy.stats.norm(0, 1).cdf((ref_mean * z - test_mean) / (ref_mean * z**2 + test_mean)**(1/2.0))
            # pn = scipy.stats.norm(0, 1).cdf((ref_mean * zn - test_mean) / (ref_mean * zn**2 + test_mean)**(1/2.0))
            #print(positive_sub)
            print(positive_sub.shape[0])
            if positive_sub.shape[0] > 0:
                data.iloc[positive_sub.index, 7] = (1 - scipy.stats.norm(0, 1).cdf((ref_mean * z - test_mean) / (ref_mean * z**2 + test_mean)**(1/2.0)))
            else:
                continue
            if negative_sub.shape[0] > 0:
                data.iloc[negative_sub.index, 7] = (scipy.stats.norm(0, 1).cdf((ref_mean * zn - test_mean) / (ref_mean * zn**2 + test_mean)**(1/2.0)))
            else:
                continue
            #data.iloc[sub.index, 6:7] = sub.iloc[:, 6:7]
    else:
        test_mean = data['test'].mean()
        ref_mean = data['ref'].mean()
        #
        global_norm = test_mean / ref_mean
        global_ratio = data['test'] / data['ref']
        data['log2'] = np.log2(data['test'] / data['ref'] / global_norm)
        positive_sub = data[data['log2'] >= 0]
        negative_sub = data[data['log2'] < 0]
        z = positive_sub['test'] / positive_sub['ref']
        zn = negative_sub['test'] / negative_sub['ref']
        # p = 1 - scipy.stats.norm(0, 1).cdf((lambda_ref * z - lambda_test) / (lambda_ref * z**2 + lambda_test)**(1/2.0))
        # pn = scipy.stats.norm(0, 1).cdf((lambda_ref * zn - lambda_test) / (lambda_ref * zn**2 + lambda_test)**(1/2.0))
        data.iloc[positive_sub.index, 7] = 1 - scipy.stats.norm(0, 1).cdf((ref_mean * z - test_mean) / (ref_mean * z ** 2 + test_mean) ** (1 / 2.0))
        data.iloc[negative_sub.index, 7] = scipy.stats.norm(0, 1).cdf((ref_mean * zn - test_mean) / (ref_mean * zn ** 2 + test_mean) ** (1 / 2.0))

    return data

#write bam subtraction results to file
def writeResults(subtractionTable):
    subtractionTable.to_csv("bam_subtraction_out_file.tab", sep="\t")
    print("BAM file subtraction is suceessfully finished")
    print("check bam_subtraction_out_file.tab file!!")

#modular function
def subtractBAMs(bam_qry, bam_ref, reference_fasta, p_value, log2, alpha, genome_size, sliding=True):

    steps = False
    if sliding:
        steps = True
    else:
        steps = False

    ref_keys, ref_start, ref_end, qry_cov, ref_cov = ([] for i in range(5))

    #This method does the following jobs:
    #check the BAM file integrity
    #create fasta dictionary
    #Calculates window size based on p-value and log2 value
    #estimate coverage based on read counts
    #returns numpy arrays for reference and query names, position, coverages
    ref_keys, ref_start, ref_end, qry_cov, ref_cov = count_coverage(bam_ref, bam_qry, reference_fasta, p_value, log2, alpha, genome_size, steps)

    #creates a pandas table with 5 columns: reference name, start, end, test, ref
    count_data = createCountData(ref_keys, ref_start, ref_end, qry_cov, ref_cov)

    #estimates log2 fold change and p-values for individual windows
    processed_data = estimateSignificance(count_data, chr_normalization=False)

    #outputs a file with processed data
    writeResults(processed_data)

def singleBam(qryBam, ref_fasta, p_value, log2, alpha, genome_size, steps=False):
    qryBamHandler = check_alignment_file(qryBam)
    fasta_map = create_fasta_map(ref_fasta)
    window_size = 50000 #get_windows(qryBamHandler, qryBamHandler, p_value, log2, alpha, genome_size)
    if steps == True:
        count_qry_map = sliding_coverage(qryBamHandler, fasta_map, window_size)
        qry_keys, qry_start, qry_end, qry_cov = extractData(count_qry_map)
    else:
        count_qry_map = non_sliding_coverage(qryBamHandler, fasta_map, window_size)
        qry_keys, qry_start, qry_end, qry_cov = extractData(count_qry_map)
    qryBamHandler.close()
    #return qry_keys, qry_start, qry_end, qry_cov
    cols = ['chromosome', 'start', 'end', 'test']
    data = pd.DataFrame({'chromosome': qry_keys, 'start': qry_start, 'end': qry_end, 'test': qry_cov}, columns=cols)
    data.to_csv("single_bam_subtraction_out_file.tab", sep="\t")
    print("BAM file subtraction is suceessfully finished")
    print("check bam_subtraction_out_file.tab file!!")

def main():
    subtractBAMs("pt_with_cnv_data.sorted.bam", "pt_with_control_data.sorted.bam", "Phaeodactylum_GCF_000150955.2_ASM15095v2_genomic.fasta", 0.01, 1.0, 1.5, 4200000, sliding=True)
    #singleBam("Danaus_unmapped_to_female_repeat_fragmented.sorted.bam", "danus_female_heterozygous_repeat_fragmented.fasta", 0.001, 1.0, 1.5, 566000000, steps=True)
    #bed_coverage("accepted_hits.bam", "d_subtract_a_with_A.bed", prefix="d_subtract_a")

if __name__ == "__main__":
    main()

##subtractBAMs("pt_with_cnv_data.sorted.bam", "pt_with_control_data.sorted.bam", "Phaeodactylum_GCF_000150955.2_ASM15095v2_genomic.fasta", 0.001, 0.6, 1.5, 27450724, sliding=True)
#https://github.com/JohnLonginotto/pybam
#pysam.sort("-o", "output.bam", "ex1.bam")
#samtools sort -o output.bam ex1.bam
#pysam.sort("-m", "1000000", "-o", "output.bam", "ex1.bam")
#tabixfile = pysam.TabixFile("example.gtf.gz")
#for gtf in tabixfile.fetch("chr1", 1000, 2000):
#print (gtf.contig, gtf.start, gtf.end, gtf.gene_id)
