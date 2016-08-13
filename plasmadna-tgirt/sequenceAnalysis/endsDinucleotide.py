#!/usr/bin/env python

from pybedtools import BedTool
from multiprocessing import Pool
import numpy as np
import pyfaidx
import sys
import os
import pandas as pd
import glob

def allDinucleotides():
    bases = ['A','C','T','G','N']
    dinucleotides = [i+j for j in bases for i in bases]
    return dinucleotides

def makeDinucleotideMatrix(header, window):
    length = window * 2
    diDict = {di : {i: 0 for i in range(length)} for di in header} 
    return diDict

def extractDinucleotide(dinucleotideCount, chrom , start, end, strand, fasta, chromLengthDict, dinucleotideHeader):
    chromLength = chromLengthDict[chrom]
    end += 1
    if start > 0 and end < chromLength:
        sequence = fasta[chrom][start:end].seq if strand == '-' else \
                fasta[chrom][start:end].reverse.complement.seq
        for i in np.arange(len(sequence)-1):
            dinucleotideCount[sequence[i:i+2]][i] += 1 
    return dinucleotideCount

def parseFragments(fragment, window, regularChrom, fasta, chromLengthDict, 
        shortDinucleotideCount, longDinucleotideCount, readEnd, dinucleotideHeader,
        smallFragments, largeFragments):
    insertSize = int(fragment.fields[4])
    chrom = fragment.chrom
    strand = fragment.strand
    if  chrom in regularChrom:
        if (readEnd == "3'" and strand == '-') or (readEnd=="5'" and strand == '+'):
            center = fragment.start
        elif (readEnd == "5'" and strand == '-') or (readEnd == "3'" and strand == '+'):
            center = fragment.end
        start = center - window
        end = center + window
        if insertSize >= smallFragments[0] and insertSize <= smallFragments[1]:
            shortDinucleotideCount = extractDinucleotide(shortDinucleotideCount, chrom, start, end, 
                    strand, fasta, chromLengthDict, dinucleotideHeader)
        elif insertSize >= largeFragments[0] and insertSize <= largeFragments[1]:
            longDinucleotideCount = extractDinucleotide(longDinucleotideCount, chrom, start, end, 
                    strand, fasta, chromLengthDict, dinucleotideHeader)
    return shortDinucleotideCount, longDinucleotideCount

def makeDF(nucleotideDict, lenType, window, fragment):
    df = pd.DataFrame.from_dict(nucleotideDict)
    df['sum'] = df.sum(axis=1)
    df['position'] = df.index - window
    df = pd.melt(df,id_vars=['position','sum'], var_name='dinucleotide', value_name = 'count')
    df['fraction'] = np.true_divide(df['count'],df['sum'])
    df['lenType'] = np.repeat(lenType, len(df))
    df.drop(['sum','count'], axis=1, inplace=True)
    return df

def parseBed(bedFile, fasta, window, chromLengthDict, regularChrom, readEnd, smallFragments, largeFragments):
    bed = BedTool(bedFile) 
    dinucleotideHeader = allDinucleotides()
    shortDinucleotideCount = makeDinucleotideMatrix(dinucleotideHeader, window)
    longDinucleotideCount = makeDinucleotideMatrix(dinucleotideHeader, window)
    for fragment in bed:
        shortDinucleotideCount, longDinucleotideCount = parseFragments(fragment, window, regularChrom, fasta, 
                chromLengthDict, shortDinucleotideCount, longDinucleotideCount, readEnd, dinucleotideHeader,
                smallFragments, largeFragments)
    iterable = zip([shortDinucleotideCount, longDinucleotideCount],
            ['short (%i-%ibp)' %(smallFragments[0], smallFragments[1]),'long (%i-%ibp)' %(largeFragments[0], largeFragments[1])],
            [smallFragments, largeFragments])
    dfs = [makeDF(nucleotideDict, lenType, window, fragment) for nucleotideDict, lenType, fragment in iterable]
    df = pd.concat(dfs)
    df['readEnd'] = np.repeat(readEnd, len(df))
    return df

def runFile(args):
    bedFile, reference, window, regularChrom, smallFragments, largeFragments = args
    sample = os.path.basename(bedFile).split('.')[0]
    print 'Start analyzing %s ' %sample
    with  pyfaidx.Fasta(reference) as fasta:
        chromLengthDict = {key:len(fasta[key]) for key in fasta.keys()}
        dfs = [parseBed(bedFile, fasta, window, chromLengthDict, 
            regularChrom, readEnd, smallFragments, largeFragments) for readEnd in ["3'","5'"]]
    df = pd.concat(dfs)
    df['samplename'] = np.repeat(sample,len(df))
    print 'Finished analyszing %s' %sample
    return df

def main():
    projectpath = '/Users/wckdouglas/plasmaDNA/results'
    bedpath = projectpath + '/bedFile'
    resultpath = projectpath  + '/dinucleotides/cleavedEnds'
    reference = '/Users/wckdouglas/plasmaDNA/reference/genome_rDNA.fa'
    window = 120
    smallFragments = (35,80)
    largeFragments = (120,180)
    if not os.path.isdir(resultpath):
        os.mkdir(resultpath)
    bedFiles = glob.glob(bedpath + '/*bed')
    regularChrom = map(str, np.arange(1,23))
    regularChrom.extend(['X','Y'])
    dfs = Pool(12).map(runFile, [(bedFile, reference, window, regularChrom, 
                                smallFragments, largeFragments) for bedFile in bedFiles]) 
    df = pd.concat(dfs)
    tablename = resultpath + '/readEndsDinucleotidesTable.tsv'
    df.to_csv(tablename, index = False, 
            index_label = False, sep = '\t')
    print 'Saved %s ' %tablename

if __name__ == '__main__':
    main()
