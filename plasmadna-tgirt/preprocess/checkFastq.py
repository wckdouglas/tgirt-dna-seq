#!/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import argparse
import gzip 
import glob

parser = argparse.ArgumentParser(description='Check splitted fastq for sequence length and quality length.')
parser.add_argument('-i','--datapath', required=True, help='Fastq path for the fastq files to be checked.')
args = parser.parse_args()
files = glob.glob(args.datapath + '/*gz')
for fastq in files:
    if 'fastq' in fastq or 'fq' in fastq:
        with gzip.open(fastq) as file_handle:
            for id, seq, qual in FastqGeneralIterator(file_handle):
                if len(seq) != len(qual):
                    print '%s: %s not matched' %(fastq, id)
        print 'Finished checking %s ' %fastq
    else:
        print '%s is not fastq file' %fastq
