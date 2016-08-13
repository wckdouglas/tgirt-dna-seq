#!/bin/env python

# This is the program for converting bed file to
# bam file and calculate window prediction scores from the selected window size 
# input: bed file for alignment, each record is a fragment
#        gtf file 
#        fai file

from pybedtools import BedTool
import numpy as np
import argparse
import pysam
import subprocess
import os
import sys
import pyBigWig as pbw
from sys import stderr, argv
from multiprocessing import Pool
from scipy.signal import savgol_filter, medfilt

def printMessage(message, sample):
    programname = os.path.basename(argv[0]).split('.')[0]
    stderr.write('[%s] %s: %s\n' %(programname,sample,message))
    return 0

def getOpt():
    parser = argparse.ArgumentParser(description='Extract coverage of TSS')
    parser.add_argument('-i','--inFile',help='input bed file (merged paired-end fragments to single fragment bed file)',required=True)
    parser.add_argument('-o','--outprefix',help='output prefix', default='out')
    parser.add_argument('-g','--genome',help='genome file (.fa.fai)',required=True)
    parser.add_argument('-c','--chromosome',help='chromosome name',required=True)
    parser.add_argument('-t','--TSSwindow', help='Window size for calculating WPS (default: 2000)', default = 2000, type=int)
    args = parser.parse_args()
    inFile = args.inFile
    outprefix = args.outprefix
    genome = args.genome
    tssWindow = args.TSSwindow
    chromosome = args.chromosome
    return inFile, outprefix, genome, tssWindow, chromosome

def makeBamFile(inFile, tempBam, samplename, genome):
    printMessage('Making %s ' %tempBam, samplename)
    command = 'bedtools bedtobam -i %s -g %s ' %(inFile, genome) +\
            '| samtools sort -T %s -O bam - ' %(tempBam.replace('.bam','')) +\
            '> %s ' %(tempBam)
    subprocess.call(command,shell=True)
    subprocess.call('samtools index %s' %tempBam, shell=True)
    printMessage('Indexed %s' %tempBam, samplename)
    return 0

def makeBam(inFile,outprefix, genome, samplename):
    '''
    bed file as input first sort then bedtobam
    index the bam file finally
    '''
    bamFolder = '/'.join(outprefix.split('/')[:-1]) +'/bamFiles'
    subprocess.call('mkdir -p %s' %bamFolder, shell=True)
    tempBam = bamFolder + '/' + outprefix.split('/')[-1] + '.sorted.bam'
    if not os.path.isfile(tempBam):
        makeBamFile(inFile, tempBam, samplename, genome)
    else:
        bamTime = os.path.getmtime(tempBam)
        bedTime = os.path.getmtime(inFile)
        if  bedTime > bamTime :
            printMessage('Old bamfile: %s... making new....' %tempBam, samplename)
            makeBamFile(inFile, tempBam, samplename, genome)
        else:
            printMessage('Used existing bamfile: %s' %tempBam, samplename)
    return tempBam


def pushWPStoArray(aln, halfWPSwindow, start, end, tssWindow, isize, wpsWindow):
    """
    for a given alignment, compute the regions that can be fully aligned and not
    e.g. [-1, -1 , -1, 1 , 1, 1, 1, 1, 1, -1, -1, -1] for a wps window -f 6 (halfwindow 3 )
    this will be added to the defined transcription start site wps score array after 
    adjusting fot the position
    """
    tssWindow = end - start
    transcriptAlnWPS = np.zeros(tssWindow) # setting the tss window as zeros wps array
    alnWPS = np.zeros(isize + wpsWindow) #adding halfwindow to both side of the alignment
    alnWPS[wpsWindow:-wpsWindow] = 1 # only half window after aln start, the alignment is fully covering the wps window 
                                     # and we added half window on the previous line
    alnWPS[alnWPS != 1 ] = -1 #making the alignment wps with ends containing windows
    alnShift = start - (aln.pos - halfWPSwindow) # the distance between alignment start and right side of the wps window: 
                                  #  + is left side of the window start, - is to the right
    if alnShift >= 0:
        wps = alnWPS[alnShift:]
        end = len(wps) if len(wps) < tssWindow else tssWindow
        transcriptAlnWPS[:end] += wps[:end]
    else:
        baseShifted = abs(alnShift)
        end = tssWindow if baseShifted + len(alnWPS) > tssWindow else baseShifted + len(alnWPS)
        alignedBases = tssWindow + alnShift
        wps = alnWPS[:alignedBases]
        transcriptAlnWPS[baseShifted:end] += wps
    return transcriptAlnWPS

def calculateWPS(bam, chrom, start, end, tssWindow, wpsWindow, halfWPSwindow, upperBound, lowerBound):
    '''
    for each gene start site region:
        for i in each position:
            fetch all alignments in the region with the window size (e.g. 120)
            calculate number of alignments fully span and partially mapped
            get wps[i] = number of full map - number of end mapped
    '''
    transcriptWps = np.zeros(end - start)
    coverage = 0
    for aln in bam.fetch(reference = chrom, start = start - halfWPSwindow, end = end+halfWPSwindow):
        isize = abs(aln.qlen)
        if upperBound > isize > lowerBound:
            transcriptWps += pushWPStoArray(aln, halfWPSwindow, start, 
                    end, tssWindow, isize, wpsWindow)
            coverage += 1
    return transcriptWps, coverage

def extractTSSaln(bam, tssWindow, wpsWindow, halfWPSwindow, upperBound, 
        lowerBound, chrom, chromSize, samplename):
    '''
    adding up wps track for all genes 
    '''
    chromArray = np.zeros(chromSize)
    for start in np.arange(halfWPSwindow, chromSize - halfWPSwindow, tssWindow):
        end = start + tssWindow if (start + tssWindow + halfWPSwindow) < chromSize else chromSize - halfWPSwindow
        wpsTSS, coverage = calculateWPS(bam, chrom, start, end, 
                    tssWindow, wpsWindow, halfWPSwindow, upperBound, lowerBound)
        chromArray[start:end] += wpsTSS
    chromArray = chromArray - medfilt(chromArray, 1001)
    chromArray = savgol_filter(chromArray, window_length = 21, polyorder=2)
    printMessage('Finished calculating WPS for chromosome %s' %(chrom), samplename)
    return chromArray

def findIntercepts(wpsArray, outFile, chromosome,samplename, lenType):
    '''
    looking for the positions that across 0 in the de-noised signal
    '''
    wpsArray = np.asarray(wpsArray)

    signs = np.sign(wpsArray)
    signs[signs==0] = -1
    start = np.where(np.diff(signs)>0)[0]  
    end = np.where(np.diff(signs)<0)[0]  
    numberOfPeak = 0
    lowPeakBound, upperPeakBound = (15, 120) if 'Long' in lenType else (50,150)
    for i, j  in zip(start, end):
        nucleosomeSize = np.abs(j - i)
        coverageScore = np.max(wpsArray[i:j])
        if nucleosomeSize > lowPeakBound and nucleosomeSize < upperPeakBound:
            numberOfPeak += 1
            peakname = '%s_peak_%i' %(chromosome,numberOfPeak)
            peakPos = (i + j)/2
            line = '\t'.join(map(str,[chromosome, i, j, peakname, peakPos ,'+', coverageScore]))
            outFile.write(line+'\n')
    printMessage('Written %i peaks to %s' %(numberOfPeak, outFile.name), samplename)
    return 0

def writeWig(chromArray, outputWig, chromosome, samplename):
    outWig =  pbw.open(outputWig,'w') 
    chrom_length = len(chromArray)
    outWig.addHeader([(chromosome,chrom_length)])
    outWig.addEntries(chromosome, range(chrom_length), values=map(np.float64,chromArray), span=1)
    outWig.close()
    printMessage('Witten %s' %outputWig, samplename)
    return 0

def runFile(arg):
    tempBam, outprefix, genome, wpsWindow, tssWindow, upperBound, lowerBound, lenType, samplename, chromosome = arg
    wpsWindow = wpsWindow + 1
    halfWPSwindow = np.divide(wpsWindow,2)
    halfTSSwindow = (tssWindow-1)/2
    outputBed = outprefix + '.'+lenType.split(' ')[0] +'.bed'
    outputWig = outputBed.replace('.bed','.bigWig')
    with pysam.Samfile(tempBam,'rb') as bam:
        chroms = np.array(bam.references)
        if chromosome not in chroms:
            sys.exit('Wrong chromosome name: %s!' %chromosome)
        chromLength = np.array(bam.lengths)
        chromSize = int(chromLength[chroms==chromosome][0])
        chromArray = extractTSSaln(bam, tssWindow, wpsWindow, halfWPSwindow, upperBound, 
                lowerBound, chromosome, chromSize, samplename)
    writeWig(chromArray, outputWig, chromosome, samplename)
    with open(outputBed,'w') as outBed:
        findIntercepts(chromArray, outBed, chromosome, samplename, lenType)
    return 0

def main(inFile, outprefix, genome, tssWindow, chromosome):
    '''
    main function for controling the work flow
    '''
    samplename = os.path.basename(inFile).split('.')[0]
    printMessage( 'Saving all result to: %s' %outprefix, samplename)
    tempBam = makeBam(inFile, outprefix, genome, samplename)
    args = [(tempBam, outprefix, genome, wpsWindow, tssWindow, upperBound, lowerBound, lenType, samplename, chromosome) \
            for upperBound, lowerBound, lenType, wpsWindow in zip([80, 180],[35, 120],['Short (35-80bp)','Long (120-180bp)'],[16,120])]
    Pool(2).map(runFile, args)
    return 0

if __name__ == '__main__':
    inFile, outprefix, genome, tssWindow, chromosome = getOpt()
    main(inFile, outprefix, genome, tssWindow, chromosome)
