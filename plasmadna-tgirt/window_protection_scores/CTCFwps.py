#!/bin/env python

# This is the program for converting bed file to
# bam file and calculate window prediction scores from the selected window size
# input: bed file for alignment, each record is a fragment
#        gtf file
#        fai file

from matplotlib import use as mpl_use
mpl_use('Agg')
from pybedtools import BedTool, set_tempdir
from pybedtools.cbedtools import Interval
import numpy as np
import argparse
import pysam
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import wpsTools
from multiprocessing import Pool
import glob
from functools import partial
from itertools import product
sns.set_style('white')

def CTCFline(feature, half_window):
    '''parsing gtf file and take all genes and protein coding records
    extract transcription start site, if strand =='-' end  is taken as start
    making a window of 1000 in both direction
    return a bed record line
    '''
    if 'CTCF' in feature.name:
        binding_site = ((feature.end) + (feature.start))/2
        start = int(binding_site) - half_window
        end = int(binding_site) + half_window
        if start > 0:
            assert (start < end), ' Wrongly parsed gene TSS'
	    CTCF = Interval(chrom = feature.chrom.replace('chr',''), 
                            start= start,
		            end = end, 
                            name = feature.name, 
                            strand = feature.strand)
	    return CTCF

def makeCTCFbed(tssFile, halfTSSwindow):
    ctcf_bed = os.path.dirname(tssFile) + '/CTCF_regions.bed'
    if not os.path.isfile(ctcf_bed):
        BedTool(tssFile)\
		.each(CTCFline,half_window = halfTSSwindow)\
		.saveas(ctcf_bed)
	print 'Saved and using: %s' %ctcf_bed
    else:
        print 'Using %s ' %ctcf_bed
    return ctcf_bed

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
            print 'Parsing %i CTCF sites' %geneCount
    return wps

def plotResult(df, figurename):
    '''
    plotting wps figure x=position, y = wps
    write table for storing data
    '''
    with sns.plotting_context('paper', font_scale=1.5):
        p = sns.FacetGrid(data = df, col = 'samplename', row = 'type', sharey=False, aspect = 1.8, size=3)
        p.map(plt.plot, 'position','wps')
        p.set_titles('{col_name} [{row_name}]')
        p.set_xlabels('Position relative to TSS')
        p.set_ylabels('Window protection score')
        p.savefig(figurename, dpi = 300)
    print 'saved %s' %(figurename)
    return 0

def runBoundaries(samplename, bedFile, bam, windowSize, ctcfBed, boundary):
    # extract boundary information
    typename = 'Short' if max(boundary) < 100 else 'Long'
    lowerBound, upperBound = boundary
    typename = typename + ' (%i-%i bp)' %(boundary[0],boundary[1])
    wpsWindow = 120 if max(boundary) > 100 else 16
    halfWPSwindow = np.divide(wpsWindow,2)
    print 'Running file: %s with insert size boundary: %s' %(bedFile,str(boundary))

    #run wpsTools function and filter signalm funally make data frame
    wps = extractTSSaln(bam, ctcfBed, windowSize, wpsWindow, halfWPSwindow, upperBound, lowerBound)
    baseline = np.mean([wps[:500], wps[-500:]])
    wps = wps - baseline
    if 'Long' in typename:
        corrected_wps = wps - pd.Series(wps)\
                                .rolling(window=200,center=True)\
                                .mean()\
                                .values
    else:
        corrected_wps = wps

    wpsDF = pd.DataFrame({'position':np.arange(windowSize)-windowSize/2,
                          'wps':wps,
                          'corrected_wps': corrected_wps}) \
        .assign(samplename = samplename) \
        .assign(type = typename)
    return wpsDF

def runFile(ctcfBed, genome, windowSize, bedFile):
    '''
    main function for controling the work flow
    '''
    boundaries = [(35,80),(120,180)]
    samplename = os.path.basename(bedFile).split('.')[0]
    tempBam = wpsTools.makeBam(bedFile, samplename, genome)
    with pysam.Samfile(tempBam,'rb') as bam:
        func = partial(runBoundaries, samplename, bedFile, bam, windowSize, ctcfBed)
        wpsDF = map(func,boundaries)
    wpsDF = pd.concat(wpsDF,axis=0)
    return wpsDF

def makedir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    return 0

def main():
    #define folders/paths
    threads = 12
    projectPath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    referencePath = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA'
    #bedPath = projectPath + '/rmdupBedFiles'
    bedPath = projectPath + '/bedFiles'
    resultPath = projectPath + '/wpsCTCF'
    makedir(resultPath)
    makedir(bedPath + '/bamFiles')
    set_tempdir(resultPath)

    #define reference files
    featuresGTF = referencePath + '/MotifFeatures.gff'
    genome = referencePath + '/genome_rDNA.fa.fai'

    #define input/output files
    bedFiles = glob.glob(bedPath + '/*.bed')
    outputprefix = resultPath + '/CTCFwps'
    figurename = outputprefix + '.pdf'
    tablename = outputprefix + '.tsv'

    windowSize = 2000 #how long to measure?

    #run files
    halfTSSwindow = (windowSize-1)/2
    ctcfBed = BedTool(makeCTCFbed(featuresGTF, halfTSSwindow))
    p = Pool(threads)
    processFileFunc = partial(runFile, ctcfBed, genome, windowSize)
    wpsDF = p.map(processFileFunc, bedFiles)
    p.close()
    p.join()
    df = pd.concat(wpsDF)
    df.to_csv(tablename, sep='\t', index=False)
    plotResult(df, figurename)
    return 0

if __name__ == '__main__':
    main()
