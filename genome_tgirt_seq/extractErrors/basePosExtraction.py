#!/usr/bin/env python

import pysam
import re
import numpy as np
import sys
from pyfaidx import Fasta
import progressbar
from multiprocessing import Pool
import glob
import os

def mismatchMatrix(seqLength):
    """
        read
    pos    A C T G N
    0      0 0 0 0 0
    1      0 0 0 0 0
    2      0 0 0 0 0
    3      0 0 0 0 0
    .      0 0 0 0 0
    .      0 0 0 0 0
    .      0 0 0 0 0
    seqlength 0 0 0 0 0
    """
    assert isinstance(seqLength,int), seqLength
    return np.zeros(seqLength * 5,dtype=np.int32).reshape(seqLength,5), {'A':0,'C':1,'T':2,'G':3,'N':4}

def cigarToSeq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: MMMMMMMMMMMMMMMMIIIIIMMMMMMMMMM
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['M','I','D'],dtype='string'))
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def MDToSeq(mdtag):
    """
    input mdTag: 31^CA31
    output: -------------------------------DDD-------------------------------

    input mdTag: 8G11^TG14T0T28T3
    output: --------G-----------DD--------------TT----------------------------T---
    """
    mdNum = np.array(re.findall('[0-9]+',mdtag),dtype=np.int16)
    mdStr = np.array(re.split('[0-9]+',mdtag),dtype='string')[1:]
    mdSeq = ''
    for s, n in zip(mdStr,mdNum):
        string = n * '-' + (len(s)-1) * 'D' if '^' in s else n * '-' + s
        mdSeq += string
    return mdSeq

def processAln(aln, misMat, tTable, qualThresh, lenCut, refSeq):
    startpos = aln.reference_start
    sequence = aln.query_alignment_sequence
    MDtag = aln.get_tag('MD')
    cigarSeq = cigarToSeq(aln.cigarstring) #(MDtags have no soft clipped)
    qual =  np.array(aln.query_alignment_qualities,dtype=np.int8)
    newSeq, newCigar, newRef, newQual = '', '','', []
    cigar_pos, seq_pos = 0, 0
    mdSeq = aln.get_reference_sequence()
    positions = aln.get_reference_positions()
    cigar_index, seq_index, ref_index = 0, 0 , 0
    while cigar_index < len(cigarSeq) and ref_index < len(mdSeq):
        c = cigarSeq[cigar_index]
        s = sequence[seq_index]
        q = qual[seq_index]
        r = mdSeq[ref_index]
        if c == 'M':
            newSeq += s
            newCigar += c
            newRef += r
            newQual.append(q)
            cigar_index += 1
            seq_index += 1
            ref_index += 1
        elif c == 'D':
            newCigar += c
            newSeq += '-'
            newRef += r
            newQual.append(0)
            cigar_index += 1
            ref_index += 1
        elif c == 'I':
            newCigar += c
            newSeq += s
            newRef += '-'
            newQual.append(q)
            cigar_index += 1
            seq_index += 1
    qlength = len(newSeq)
    assert len(newRef) == qlength, '\n' + MDtag +\
                            '\n'+ aln.cigarstring+\
                            '\n'+sequence + \
                            '\nref  ' + newRef + \
                            '\nseq: ' + newSeq + \
                            '\ncig: ' + newCigar +\
                            '\nqul: ' + ''.join(map(lambda q: chr(q+33), newQual))
    newSeq, newRef, newCigar = map(lambda x: np.array(list(x)),[newSeq,newRef,newCigar])
    boolean = (newRef != '-') & (newSeq != '-')
    newSeq, newRef, newCigar = map(lambda x: ''.join(map(str,x[boolean])),[newSeq,newRef,newCigar])
    assert len(newRef) == len(positions), 'ref seq len: %i \npositions len: %i' %(len(newRef),len(positions)) +\
        '\nref: %s' %(newRef) +\
        '\nred: %s' %(newSeq) +\
        '\ncig: %s' %(newCigar) +\
        '\n %s' %(aln.cigarstring)
    pos = np.arange(qlength)
    upper = qlength - lenCut
    for b, m, q, p, rp in zip(newSeq, newRef, newQual, pos, positions):
        if q > qualThresh and c == 'M' and upper >= p >= lenCut:
            misMat[rp,tTable[b]] += 1
    return 0

def generate_intervals(refsize, window):
    sample_num = refsize / float(window)
    intervals = np.linspace(0,refsize, int(sample_num)) if sample_num >= 2 else [0,refsize]
    intervals = np.unique(np.array(intervals,dtype=np.int64))
    intervals = [(intervals[i],intervals[i+1]) for i in np.arange(len(intervals)-1)]
    return intervals

def parseRegion(refID, refsize, bam, misMat, tTable, qualThresh, lenCut, refSeq, out_file, window):
    count = 0
    intervals = generate_intervals(refsize, window)
    job_length = len(intervals)
    bar =  progressbar.ProgressBar(maxval=job_length,
	        widgets=[progressbar.Bar('=', '[', ']'),
		    ' ', progressbar.Percentage()])
    sys.stderr.write('For chromosome: %s:\n' %refID)
    bar.start()
    status = 0
    for i, interval in enumerate(intervals):
        start, end = interval
        for aln in bam.fetch(refID, start, end) :
            if not aln.is_supplementary and not aln.is_secondary and aln.is_proper_pair:
                processAln(aln, misMat, tTable, qualThresh, lenCut, refSeq)
            count += 1
            if count % 100000 == 0:
                sys.stderr.write('Processed %i alignments in %s\n' %(count, refID))
        bar.update(status + i)
    bar.finish()
    rownum = 0
    for rownum, row in enumerate(misMat):
        cov =  sum(row)
        if cov > 0:
            ref_base = str(refSeq[rownum])
            file_line =  '\t'.join([refID, str(rownum), ref_base, str(cov), '\t'.join(map(str,row))])
            out_file.write(file_line + '\n')
    return 0


def parse_bam(inBam, refFa, qualThresh, lenCut, out_path, window):
    sys.stderr.write('Starting %s\n' %(inBam))
    indexing =  pysam.index(inBam)
    samplename = os.path.basename(inBam).split('.')[0]
    with pysam.Samfile(inBam,'rb') as bam, Fasta(refFa) as fa:
        refnames = bam.references
        refsizes = bam.lengths
        for refID, refsize in zip(refnames, refsizes):
            out_file_name = out_path + '/' + samplename + '_bases_' + refID + '.tsv'
            with open(out_file_name,'w') as out_file:
                out_file.write('refID\trefPos\trefBase\tcoverage\tA\tC\tT\tG\tN\n')
                misMat, tTable = mismatchMatrix(refsize)
                refSeq = fa[refID]
                parseRegion(refID, refsize, bam, misMat, tTable, qualThresh, lenCut, refSeq, out_file, window)
            sys.stderr.write('Parsed %s in %s\n' %(refID, inBam))
    sys.stderr.write('Finished parsing %s\n' %(inBam))
    return 0

def main():
    project_path = '/scratch/02727/cdw2854/jurkatCells'
    out_path = project_path + '/pileup_bases'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    refFa = '/corral-repl/utexas/2013lambowitz/Ref/hg19/Sequence/WholeGenomeFasta/genome.fa'
    qualThresh = 33
    lenCut = 3
    window = 10000
    if len(sys.argv) != 2:
        sys.exit('usage: python %s <bamfile>' %(sys.argv[0]))
    inBam = sys.argv[1]
    parse_bam(inBam, refFa, qualThresh, lenCut, out_path, window)

if __name__=='__main__':
    main()
