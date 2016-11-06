#!/usr/bin/env python

from __future__ import print_function
from pyfaidx import Fasta
import numpy as np
import os
import string
import sys
if len(sys.argv) != 2:
    sys.exit('usage: python %s <ref_fasta>' %sys.argv[0])

complement = string.maketrans('ACTGNactgn','TGACNTGACN')

def reverse_complement(seq):
    return seq.translate(complement)[::-1]


reference = Fasta(sys.argv[1])
for line in sys.stdin:
    chrom, start, end, name, length, strand = line.rstrip().split('\t')
    start, end, length =  map(int, [start, end, length]) 
    sequence = str(reference.get_seq(chrom, start + 1, end + 1))
    if length > 150:
        read1 = sequence[:150]
        read2 = reverse_complement(sequence[-150:])
        qual = 'I' * 150
    else:
        read1 = sequence
        read2 = reverse_complement(sequence)
        qual = 'I' * (length + 1)
    line1 = '@%s/1\n%s\n+\n%s\n' %(name, read1, qual)  
    line2 = '@%s/2\n%s\n+\n%s' %(name, read2, qual)
    line = line1 + line2
    print(line, file=sys.stdout)
