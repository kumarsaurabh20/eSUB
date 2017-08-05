#!/usr/bin/env python2
from __future__ import print_function
import sys
import gzip

# Splitting the Reads
# -------------------
# This is optional, but it provides for a way to fully utilize a computing
# cluster. If we have a cluster of N nodes and the genome was split into Y pieces,
# we would want to split the reads into N/Y pieces. In our example, say N=4 and
# Y=2, so we split the reads into 2 pieces, qr-1of2.fa and qr-2of2.fa.
# To split reads, use the splitreads.py script.

def handle_gz_or_ascii(file):
    fd = None
    try:
        fd = gzip.GzipFile(file)
    except IOError:
        fd = open(file, "r")
    try:
        line = fd.readline()
    except IOError:
        fd.close()
        fd = open(file, "r")
    fd.seek(0)
    return (fd)

if sys.argv[1] == '--paired' or sys.argv[1] == '-p':
    mode = 'paired'
    mod = 1
else:
    mode = 'unpaired'
    mod = 0

READS_PER_FILE = int(sys.argv[1 + mod])

if mode == 'paired':
    READS_PER_FILE *= 2

out = None
readsdone = 0

suffix = "fasta"
if sys.argv[2].lower().strip().endswith(".fasta"):
    suffix = "fasta"
elif sys.argv[2].lower().strip().endswith(".fsa"):
    suffix = "fsa"
elif sys.argv[2].lower().strip().endswith(".fa"):
    suffix = "fa"

fd = handle_gz_or_ascii(sys.argv[2 + mod])
for line in fd:
    if line.startswith("#"):
        continue
    if line.startswith(">") and (readsdone % READS_PER_FILE) == 0:
        if out != None:
            out.close()

        fname = "%u_to_%u.%s" % (readsdone, readsdone + READS_PER_FILE - 1, suffix)

        out = open(fname, "w")
    out.write(line)

    if line.startswith(">"):
        readsdone = readsdone + 1

if out != None:
    out.close()
