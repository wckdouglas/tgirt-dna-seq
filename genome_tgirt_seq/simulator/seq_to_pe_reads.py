#!/usr/bin/env python

from __future__ import print_function
from pyfaidx import Fasta
import numpy as np
import os
import string
import sys
if len(sys.argv) != 2:
    sys.exit('usage: python %s <read_length>' %sys.argv[0])

complement = string.maketrans('ACTGNactgn','TGACNTGACN')

def reverse_complement(seq):
    return seq.translate(complement)[::-1]


read_length = int(sys.argv[1])

' Input must be tab deliminated: name   seq'

for line in sys.stdin:
    name, sequence = line.rstrip().split('\t')
    length = len(sequence)
    if length > read_length:
        read1 = sequence[:read_length]
        read2 = reverse_complement(sequence[-read_length:])
        qual = 'I' * read_length
    else:
        read1 = sequence
        read2 = reverse_complement(sequence)
        qual = 'I' * length

    line1 = '@%s/1\n%s\n+\n%s\n' %(name, read1, qual)
    line2 = '@%s/2\n%s\n+\n%s' %(name, read2, qual)
    line = line1 + line2
    print(line, file=sys.stdout)
