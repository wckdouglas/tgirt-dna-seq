#!/bin/env python

import pysam #conda: pysam 0.6
import sys
import numpy as np
import re
import os
import argparse
from multiprocessing import Pool, Manager
import progressbar
programnam = sys.argv[0]

def getOption():
    """
    Get all option fomr command lines
    """
    parser = argparse.ArgumentParser(description='Pile up bam file and get mapped bases (excluded clipped bases and low quality base).')
    parser.add_argument('-i','--bamfile',required=True,
            help='position sorted bam file')
    parser.add_argument('-q','--qualThresh',type=int, default = 33,
            help='Base calling quality threshold (default: 33)')
    parser.add_argument('-r','--refFasta', required=True,
            help='reference fasta')
    parser.add_argument('-d','--depth', type=int, 
            default = 8000, help='Maximum depth (default: 8000)')
    parser.add_argument('-p','--threads', type=int, 
            default = 1, help='Threads to used (default: 1)')
    parser.add_argument('-s','--skipBases', type=int,  
            default = 3, help='Bases at the end positions that will not be evaluated (default: 3)')
    parser.add_argument('-o','--output', type = argparse.FileType('w'),  
            default = sys.stdout, help='Output file')
    args = parser.parse_args()
    return args
    
def printLine(chrom_name,referenceBase, bases, count , position, coverage, outFile):
    """
    for each column position (pileup column), output the 
    1. position on the chromosome
    2. reference base on the ref-fasta
    3. count of A, T, C, G
    """
    regularBase = np.array(['A','T','C','G'],dtype='string')
    coverage = np.sum(count[np.in1d(bases,regularBase)]) #only high qual base were counted in coverage
    baseDict = {base: count[bases==base][0] if base in bases else 0 for base in regularBase}
    line = '\t'.join([chrom_name,str(position),str(position+1),referenceBase,str(coverage)])
    for base in regularBase:
        line += '\t' + str(baseDict[base])
    if coverage > 0:
        outFile.write(line+'\n')
    return 0

def cigarToSeq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: SMMMMMMMMMMMMMMMMIIIIIMMMMMMMMMMSSSS
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['S','M','I'],dtype='string'))
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def extractBase(args):    
    """
    for each alignment that mapped to the position,
    extract the base and the cigar annotation
    if the base is a Mapped position ('M') and have quality higher than the 
    given threshold, return the base
    """
    sequence, quality, cigar, matchPos, qualThresh, skipBases = args
    if matchPos > (skipBases-1) and matchPos < (len(sequence)-skipBases) and not np.any(np.in1d(['I','D'],list(cigar))):
        cigarSeq = cigarToSeq(cigar)
        assert len(cigarSeq) == len(sequence), '\n%s\n%s' %(cigarSeq,sequence)
        qual = quality[matchPos]
        base = sequence[matchPos]
        if (ord(qual) - 33) > qualThresh and cigarSeq[matchPos] == 'M':
	    return base

def analyzePosition(chrom_name ,pileupColumn, refBase, threads, position, qualThresh, skipBases, outFile):
    """
    for each pileup position, extracted all alignments using pysam
    processing each alignment is computationally heavy.
    used multiprocessing in here
    """
    cov = pileupColumn.n
    pool = Pool(threads)
    result = pool.map_async(extractBase, [(aln.alignment.seq, aln.alignment.qual, 
                                                    aln.alignment.cigarstring, aln.query_position,
                                                    qualThresh, skipBases) for aln in pileupColumn.pileups \
                                    if (not aln.alignment.is_secondary and aln.indel==0)])
    pool.close()
    pool.join()
    posbases = filter(None, result.get())
    bases, count = np.unique(posbases,return_counts=True)
    printLine(chrom_name, refBase, bases, count, position, cov, outFile) 
    return 0

def printHeader(outfile):
    outfile.write('chrom\tstart\tend\tref\tcov\tA\tT\tC\tG\n')
    return 0

def parse_chrom(chrom_name, chrom_len, bam, refFasta, depth, programnam, threads, qualThresh, skipBases, outFile):
    chrom_name = int(chrom_name) if chrom_name.isdigit() else chrom_name
    for pileupColumn in bam.pileup(chrom_name,fastafile=refFasta,max_depth=depth):
        position = pileupColumn.pos
        refBase = refFasta.fetch(chrom_name, position, position+1)
        analyzePosition(chrom_name, pileupColumn, refBase, threads, position, qualThresh, skipBases, outFile)
    return 0

def main(args): 
    """
    Using pysam to pileup bam file
    1. extract aligned reads at each position
    2. extract the base on each mapped aligned reads on the position (computationally heavy)
    3. write out base count lines.
    """
    bamfile = args.bamfile
    qualThresh = args.qualThresh
    ref = args.refFasta
    depth =  args.depth
    threads =  args.threads
    skipBases = args.skipBases
    outFile = args.output

    # start program
    refFasta = pysam.Fastafile(ref) 

    #index bam
    index = bamfile + '.bai'
    if os.path.exists(index):
        os.remove(index)
        sys.stderr.write('[%s] Removed original index: %s \n' %(programnam,index))
    sys.stderr.write('[%s] Indexing %s \n' %(programnam,bamfile))
    pysam.index(bamfile)
    printHeader(outFile)
    with pysam.Samfile(bamfile,'rb') as bam:
        sys.stderr.write('[%s] Pileup BAM file: %s \n' %(programnam,bamfile))
        sys.stderr.write('[%s] Using threads:   %i \n' %(programnam,threads))
        chrom_names = bam.references
        chrom_lens = bam.lengths
	bar = progressbar.ProgressBar(maxval=len(chrom_names),
		widgets=[progressbar.Bar('=', '[', ']'), ' ', 
		    progressbar.Percentage()]) 
	bar.start()
	status = 0
	bar.update(status)
        for chrom_name, chrom_len in zip(chrom_names, chrom_lens):
            parse_chrom(chrom_name, chrom_len, bam, refFasta, depth, programnam, threads, qualThresh, skipBases, outFile)
	    status += 1
	    bar.update(status)
	bar.finished()
    sys.stderr.write('[%s] Finished extracting %s.\n' %(programnam,bamfile))
    return 0
        
if __name__ == '__main__':
    args = getOption()
    main(args)
