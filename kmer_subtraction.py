#!/usr/bin/env python2
from __future__ import print_function
from Bio import SeqIO
from Bio.Seq import reverse_complement
from multiprocessing import Pool
import datetime, sys, os, glob, resource

##database storage issue
#https://stackoverflow.com/questions/2376846/which-key-value-store-is-the-most-promising-stable
##https://docs.python.org/2/library/shelve.html
##http://www.zodb.org/en/latest/tutorial.html
##http://api.mongodb.com/python/current/tutorial.html
##https://www.mongodb.com/blog/post/getting-started-with-python-and-mongodb
##http://www.bogotobogo.com/python/MongoDB_PyMongo/python_MongoDB_pyMongo_tutorial_installing.php
##https://realpython.com/blog/python/introduction-to-mongodb-and-python/
##https://marcobonzanini.com/2015/09/07/getting-started-with-mongodb-and-python/
##http://codehandbook.org/pymongo-tutorial-crud-operation-mongodb/

##storage in mongoDB
##https://stackoverflow.com/questions/4327723/mongodb-limit-storage-size

#https://stackoverflow.com/questions/552744/how-do-i-profile-memory-usage-in-python
# def using(point=""):
#     usage=resource.getrusage(resource.RUSAGE_SELF)
#     return '''%s: usertime=%s systime=%s mem=%s mb
#            '''%(point,usage[0],usage[1],
#                 (usage[2]*resource.getpagesize())/1000000.0 )

#Update the code
#First sort the main fasta file
#make sure it has equal number of sequences
#divide the file in to smaller chunks (chunk size is dependent on file size)
#perform the setoperations on each chunks
#final sets should be merged and perform final set operations
#write the headers


#https://stackoverflow.com/questions/28859095/most-efficient-method-to-check-if-dictionary-key-exists-and-process-its-value-if
def convertBinaryInt(binary):
	return int('1' + str(binary), 2)

def convertIntBinary(num):
	return format(int(num), 'b')[1:]

#check if the file is zipped or not
#check if the input file is fastq or fasta
#if fastq then convert from fastq to fasta
#return some basic stats on sequence files
def filecheck():
	print("filecheck")

def batch_iterator(iterator, batch_size):
	entry = True
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.next()
			except StopIteration:
				entry = None
			if entry is None:
				break
			batch.append(entry)
		if batch:
			yield batch

def splitData(infile, prefix, batchsize=10000):
	filenames = []
	record_iter = SeqIO.parse(open(infile), "fasta")
	for i, batch in enumerate(batch_iterator(record_iter, batchsize)):
		filename = "%s_group_%i.fasta"%(prefix, i + 1)
		filenames.append(filename)
		with open(filename, "w") as handle:
			count = SeqIO.write(batch, handle, "fasta")
		print("Wrote %i records to %s"%(count, filename))
	return filenames

def codeSequence(fasta_file):
	kmer_size=15
	currentTime = getDateTimeString()
	fastafilename, fastafileextension = os.path.splitext(fasta_file)
	kmerBox = set()
	seqCounter = 0
	print("Reading file %s"%(fastafilename))
	for e in SeqIO.parse(fasta_file, "fasta"):
		seqdb1 = []
		sequence = e.seq
		header = e.id
		seqCounter += 1
		counter = 0
		for i in range(0, len(sequence) - kmer_size + 1):
			seq_15 = sequence[i:i + kmer_size]
			seq_15.strip()  # https://stackoverflow.com/questions/9347419/python-strip-with-n
			seq_15.upper()
			# print(seq_15)
			seq_15_rc = reverse_complement(seq_15)
			seq15bin = ''
			if str(seq_15) < str(seq_15_rc):
				seq15bin = seq_15
			else:
				seq15bin = seq_15_rc
			temp = ''
			counter += 1
			for e in str(seq15bin):
				if e == "A":
					temp += '00'
				elif e == "T":
					temp += '11'
				elif e == "G":
					temp += '10'
				elif e == "C":
					temp += '01'
				else:
					# print("Unknown character found!!")
					continue
			seq15dec = convertBinaryInt(temp)
			seqdb1.append(seq15dec)
			#kmerMap[header] = seqdb1
			kmerBox.update(seqdb1)
			# print("Read# %d"%counter1)
			print("binned %i kmer from sequence %s" %(counter, header), end="\r")
			sys.stdout.flush()
		print("Processed %i sequences from fasta"%(seqCounter), end="\r")
		sys.stdout.flush()
	kmerboxoutfile = "set" + "_" + fastafilename + "_" + currentTime + "." + "list"
	writeFile(kmerBox, kmerboxoutfile, "kmer-temps")
	#writeFile(kmerMap, "Kmer_header_database.tab", "kmer-header-temps")

def getDateTimeString():
	return datetime.datetime.now().strftime("%Y%m%d%H%M%S")

def setOperations(kmerset1, kmerset2):
	subtractedSet = (kmerset1 - kmerset2)
	commonSet = (kmerset1 & kmerset2)
	exclusiveSet = (kmerset1 ^ kmerset2)
	return subtractedSet, commonSet, exclusiveSet

def kmerHashFilter(kmer_map, subtracted_set, proportion=0.5):
	distinct_header_list=[]
	headerCount = 0
	for k,v in kmer_map.items():
		headerCount += 1
		counter = 0
		for e in v:
			if e in subtracted_set:
				counter += 1
			else:
				continue
			if float(counter/len(v)) >= float(proportion):
				distinct_header_list.append(k)
			else:
				continue
	return set(distinct_header_list)

def decodeKmers(dataset):
	diffContainer = []
	print("Decoding from binary to nucleotides!!")
	for seq in dataset:
		dummy = convertIntBinary(seq)
		oriseq = ''
		for e in range(0, len(dummy), 2):
			bits = dummy[e:e + 2]
			if bits == "00":
				oriseq += 'A'
			elif bits == "11":
				oriseq += 'T'
			elif bits == "10":
				oriseq += 'G'
			elif bits == "01":
				oriseq += 'C'
			else:
				print("Unknown character found!!")
				sys.exit(0)
		diffContainer.append(oriseq)
	return diffContainer

def writeFile(data, filename, job):
	out = open(filename, "a")
	print("Writing %s file.." % (job))
	try:
		if isinstance(data, (list, set, tuple)):
			for each in dataArray:
				out.write("%s\n" % (each))
		elif isinstance(data, dict):
			for k,v in data.items():
				out.write("%s\t%s\n"%(k, ' '.join(v)))
		else:
			print("Could not identify data type!!")
	except TypeError:
		print("Oops!  That was no valid data type.  Try again...")
	
	out.close()
	print("Finished Writing %s" % (job))

def checkMemory(set1, set2, set3, kmerMap):
	#####################Troubleshoot###################
	print("#####################")
	print("Set1 in kB for 1 read")
	print(sys.getsizeof(set1))
	print("Set2 in kB for 1 read")
	print(sys.getsizeof(set2))
	print("Set3 in kB for 1 read")
	print(sys.getsizeof(set3))
	print("kmer map size in kB for 1 read")
	print(sys.getsizeof(kmerMap))
	print("######################")
	print(sys.getsizeof(100111010101010101010110110110))
	print(sys.getsizeof("TGCGATAGACAGTAGACAGATGACAGATGA"))
	#####################Troubleshoot###################

def main():
	header_map = {}
	kmer_map = {}
	currentTime = getDateTimeString()
	filenames = splitData("Phaeodactylum_GCF_000150955.2_ASM15095v2_genomic_data_1.fa", "reference", 500000)
	p = Pool(4)
	p.map(codeSequence, filenames)
	#writeFile(kmer_map, "kmer_headers_databank.tab", "kmer-headers-temps")


if __name__ == "__main__":
	main()
