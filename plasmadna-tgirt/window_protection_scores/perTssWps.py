#!/bin/env python

# This is the program for converting bed file to
# bam file and calculate window prediction scores from the selected window size 
# input: bed file for alignment, each record is a fragment
#        gtf file 
#        fai file

import matplotlib
matplotlib.use('Agg')
from pybedtools import BedTool
import numpy as np
import argparse
import pysam
import subprocess
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
from scipy.signal import argrelextrema
import wpsTools
sns.set_style('white')

def gtfToDict(gtfFields):
    d = {}
    for item in gtfFields[-1].strip(';').split(';'):
        item = item.strip(' ')
        itemList = item.replace('"','').split(' ')
        d[itemList[0].strip(' ')] = itemList[1].strip(' ')
    return d

def gtf2Bedline(gtfLine, halfTSSwindow): 
    '''parsing gtf file and take all genes and protein coding records
    extract transcription start site, if strand =='-' end  is taken as start
    making a window of 1000 in both direction
    return a bed record line
    '''
    gtfFields = gtfLine.fields
    if gtfFields[2] == 'gene':
        if gtfLine.strand == '+':
            start = int(gtfLine.start) - halfTSSwindow
            end = int(gtfLine.start) + halfTSSwindow
        else:
            start = int(gtfLine.end) - halfTSSwindow
            end = int(gtfLine.end) + halfTSSwindow
        if start > 0:
            assert (start < end), ' Wrongly parsed gene TSS'
            infoDict = gtfToDict(gtfFields)
            name =infoDict['gene_name'] if 'gene_name' in infoDict.keys() else infoDict['gene_id']
            bedline = map(str,[gtfLine.chrom, 
                        start, end,
                        infoDict['gene_id'], 0,
                        gtfLine.strand, 
                        name,
                        infoDict['gene_biotype']])
            return '\t'.join(bedline) + '\n'

def makeTSSbed(tssFile, halfTSSwindow):
    '''
    feed in gtf file to gtf2Bedline
    optain string of bed record lines
    and generate a bed file
    '''
    geneFile = os.path.dirname(tssFile) + '/tss.all.bed'
    if not os.path.isfile(geneFile):
        gtfFile = BedTool(tssFile)
        tssLines = [gtf2Bedline(gtfLine, halfTSSwindow) for gtfLine in gtfFile]
        tssLines = filter(None,tssLines)
        geneCount = len(tssLines)
        tssLines = ''.join(tssLines)
        BedTool(tssLines,from_string=True).saveas(geneFile)
        print 'saved %s with %i genes' %(geneFile,geneCount)
    else:
        geneCount = sum(1 for line in open(geneFile,'ru'))
        print 'Using %s with %i genes' %(geneFile, geneCount)
    return geneFile, geneCount


def AutoCorrelation(x):
    x = np.asarray(x)
    y = x-x.mean()
    result = np.correlate(y, y, mode='full')
    result = result[len(result)//2:]
    result /= result[0]
    return result 

def extractPeriodicity(wps):
    corLags = AutoCorrelation(wps)
    maximas = argrelextrema(corLags, np.greater)[0]
    return 0 if maximas.size == 0 else maximas[0] 

def extractTSSaln(bam, tssBed, tssWindow, wpsWindow, halfWPSwindow, geneCount, outprefix, upperBound, lowerBound):
    '''
    adding up wps track for all genes 
    '''
    geneName = np.zeros(geneCount, dtype=object)
    periodicity = np.zeros(geneCount,dtype=np.int64)
    geneID = np.zeros(geneCount, dtype=object)
    geneType = np.zeros(geneCount,dtype=object)
    i = 0
    for tss in tssBed:
        geneName[i] = tss.fields[-2]
        geneID[i] = tss.name
        geneType[i] = tss.fields[-1]
        wps, coverage = wpsTools.calculateWPS(bam, tss, tssWindow, wpsWindow, halfWPSwindow, upperBound, lowerBound)
        if np.sum(wps) != 0 and coverage > 100:
            periodicity[i] = extractPeriodicity(wps)
        else:
            periodicity[i] = -1
        if i % 10000 == 0:
            print 'Parsing %i genes' %i
        i += 1
    dfMatrix = np.append(geneID, [geneName, geneType,periodicity]).reshape(4,geneCount).transpose()
    df = pd.DataFrame(dfMatrix,columns=['geneID','geneName','geneType','periodicity'])
    tablename = outprefix + '.tsv'
    df.to_csv(tablename,sep='\t', index = False)
    print 'saved %s' %(tablename)
    return 0

def main(inFile, outprefix, tssFile, genome, wpsWindow, tssWindow, upperBound, lowerBound):
    '''
    main function for controling the work flow
    '''
    print 'Saving all result to: %s' %outprefix
    halfWPSwindow = np.divide(wpsWindow,2)
    halfTSSwindow = (tssWindow-1)/2
    geneFile, geneCount = makeTSSbed(tssFile, halfTSSwindow)
    tssBed = BedTool(geneFile)
    tempBam = wpsTools.makeBam(inFile, outprefix, genome)
    with pysam.Samfile(tempBam,'rb') as bam:
        extractTSSaln(bam, tssBed, tssWindow, wpsWindow, halfWPSwindow, geneCount, outprefix, upperBound, lowerBound)
    return 0

if __name__ == '__main__':
    inFile, outprefix, tssFile, genome, wpsWindow, tssWindow, upperBound, lowerBound = getOpt()
    main(inFile, outprefix, tssFile, genome, wpsWindow, tssWindow, upperBound, lowerBound)
