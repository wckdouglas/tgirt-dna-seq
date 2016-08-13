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
import pylab as plt
import seaborn as sns
import wpsTools
import glob
from multiprocessing import Pool
from scipy.signal import medfilt
from scipy.signal import savgol_filter
sns.set_style('white')

def gtfToDict(gtfFields):
    infoDict = {}
    fields = gtfFields[8]
    for info in fields.split(';')[:-1]:
        geneInfo = info.replace('"','').split(' ')
        infoDict[geneInfo[0]] = geneInfo[1]
    return infoDict

def gtf2Bedline(gtfLine, halfTSSwindow, chromosomes): 
    '''parsing gtf file and take all genes and protein coding records
    extract transcription start site, if strand =='-' end  is taken as start
    making a window of 1000 in both direction
    return a bed record line
    '''
    gtfFields = gtfLine.fields
    if gtfFields[2] == 'gene' and gtfFields[1] == 'protein_coding' and gtfLine.chrom in chromosomes:
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
                        name])
            return '\t'.join(bedline) + '\n'

def makeTSSbed(genesGTF, halfTSSwindow):
    '''
    feed in gtf file to gtf2Bedline
    optain string of bed record lines
    and generate a bed file
    '''
    tssBed = os.path.dirname(genesGTF) + '/tss.bed'
    if not os.path.isfile(tssBed):
        chromosomes = map(str,np.arange(1,23))
        chromosomes.extend(['X','Y'])
        gtfFile = BedTool(genesGTF)
        tssLines = [gtf2Bedline(gtfLine, halfTSSwindow, chromosomes) for gtfLine in gtfFile]
        tssLines = filter(None,tssLines)
        tssLines = ''.join(tssLines)
        BedTool(tssLines,from_string=True).saveas(tssBed)
        print 'saved %s' %tssBed
    else:
        print 'Using %s' %tssBed
    return tssBed

def extractTSSaln(bam, tssBed, windowSize, wpsWindow, halfWPSwindow, upperBound, lowerBound):
    '''
    adding up wps track for all genes 
    '''
    geneCount = 0
    wps = np.zeros(windowSize)
    for tss in tssBed:
        geneCount += 1
        wpsTSS, coverage = wpsTools.calculateWPS(bam, tss, windowSize, wpsWindow, halfWPSwindow, upperBound, lowerBound)
        wps += wpsTSS
        if geneCount % 5000 == 0:
            print 'Parsing %i genes' %geneCount
    return wps

def plotResult(df, figurename):
    '''
    plotting wps figure x=position, y = wps
    write table for storing data
    '''
    with sns.plotting_context('paper', font_scale=1.5):
        p = sns.FacetGrid(data = df, row = 'samplename', col = 'type', sharey=False, aspect = 1.8, size=3) 
        p.map(plt.plot, 'position','wps')
        p.set_titles('{row_name} [{col_name}]')
        p.set_xlabels('Position relative to TSS')
        p.set_ylabels('Window protection score')
        p.savefig(figurename, dpi = 300)
    print 'saved %s' %(figurename)
    return 0
    
def runFile(args):
    '''
    main function for controling the work flow
    '''
    bedFile, boundary, resultPath, tssBed, genome, windowSize = args
    samplename = os.path.basename(bedFile).split('.')[0]
    typename = 'Short' if max(boundary) < 100 else 'Long' 
    typename = typename + ': (%i-%i bp)' %(boundary[0],boundary[1]) 
    wpsWindow = 120 if max(boundary) > 100 else 16
    print 'Running file: %s with insert size boundary: %s' %(bedFile,str(boundary))
    halfWPSwindow = np.divide(wpsWindow,2)
    tempBam = wpsTools.makeBam(bedFile, samplename, genome)
    lowerBound, upperBound = boundary
    with pysam.Samfile(tempBam,'rb') as bam:
        wps = extractTSSaln(bam, tssBed, windowSize, wpsWindow, halfWPSwindow, upperBound, lowerBound)
    wps = wps - medfilt(wps,201)
    wpsDF = np.array([np.arange(windowSize)-windowSize/2, wps]).transpose()
    wpsDF = pd.DataFrame(wpsDF,columns=['position','wps'])
    wpsDF['samplename'] = np.repeat(samplename, len(wpsDF))
    wpsDF['type'] = np.repeat(typename, len(wpsDF))
    return wpsDF

def makedir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    return 0

def main():
    #define folders/paths
    projectPath = '/scratch/cdw2854/plasmaDNA'
    referencePath = '/scratch/cdw2854/plasmaDNA/reference'
    bedPath = projectPath + '/bedFiles'
    resultPath = projectPath + '/wpsTSS'
    makedir(resultPath)

    #define reference files
    genesGTF = referencePath + '/genes.gtf'
    genome = referencePath + '/genome_rDNA.fa.fai'

    #define input/output files
    bedFiles = glob.glob(bedPath + '/*bed')
    outputprefix = resultPath + '/tssWPS'
    figurename = outputprefix + '.pdf'
    tablename = outputprefix + '.tsv'

    boundaries = [(35,80),(120,180)]
    windowSize = 2000 #how long to measure?
    #run files
    halfTSSwindow = (windowSize-1)/2
    tssBed = BedTool(makeTSSbed(genesGTF, halfTSSwindow))
    args = [(bedFile, boundary, resultPath, tssBed, genome, windowSize) for bedFile in bedFiles for boundary in boundaries]
    wpsDF = Pool(24).map(runFile, args)
    df = pd.concat(wpsDF)
    df.to_csv(tablename, sep='\t', index=False)
    plotResult(df, figurename)
    return 0

if __name__ == '__main__':
    main()
