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
