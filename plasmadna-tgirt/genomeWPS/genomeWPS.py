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
from scipy.signal import savgol_filter
import pandas as pd
from functools import partial
from itertools import izip

def make_folder(folder):
    if not os.path.isdir(folder):
        os.mkdir(folder)
    return 0


def printMessage(message, sample):
    programname = os.path.basename(argv[0]).split('.')[0]
    stderr.write('[%s] %s: %s\n' %(programname,sample,message))
    return 0

def getOpt():
    parser = argparse.ArgumentParser(description='Extract coverage of TSS')
    parser.add_argument('-i','--inFile',
            help='input bed file (merged paired-end fragments to single fragment bed file)',
            required=True)
    parser.add_argument('-o','--outprefix',help='output prefix', default='out')
    parser.add_argument('-g','--genome',help='genome file (.fa.fai)',required=True)
    parser.add_argument('-c','--chromosome',help='chromosome name',required=True)
    parser.add_argument('-t','--TSSwindow', help='Window size for calculating WPS (default: 2000)',
            default = 2000, type=int)
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
    for aln in bam.fetch(reference = chrom,
                         start = start - halfWPSwindow,
                         end = end+halfWPSwindow):
        isize = abs(aln.qlen)
        if upperBound > isize > lowerBound:
            transcriptWps += pushWPStoArray(aln, halfWPSwindow, start,
                    end, tssWindow, isize, wpsWindow)
            coverage += 1
    return transcriptWps, coverage


def adjust_median(wps_array):
    rolling_median = pd.Series(wps_array)\
            .rolling(window = 1000)\
            .median()
    adjusted_wps = np.nan_to_num(wps_array - rolling_median)
    return adjusted_wps


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
    chromArray = savgol_filter(adjust_median(chromArray), window_length = 21, polyorder=2)
    return chromArray


def find_peak_region(wpsArray):
    wpsArray = np.asarray(wpsArray)
    signs = np.sign(wpsArray)
    signs[signs==0] = -1
    start = np.where(np.diff(signs)>0)[0]
    end = np.where(np.diff(signs)<0)[0]
    return start, end


def merge_peaks(peak_start, peak_end):
    new_start = []
    new_end = []
    tolerance_unprotected = 5
    i = 0
    while i < len(peak_start)-2:
        new_start.append(peak_start[i])
        j = i
        while peak_start[j+1] - peak_end[j] <= tolerance_unprotected:
            j += 1
        new_end.append(peak_end[j])
        j += 1
        i = j
    new_start.append(peak_start[i])
    new_end.append(peak_end[i])
    return np.array(new_start), np.array(new_end)

def maxSubArray(ls):
    '''
    #https://gist.github.com/alabid/3734606
    '''
    if len(ls) == 0:
       raise Exception("Array empty") # should be non-empty

    runSum = maxSum = ls[0]
    i = 0
    start = finish = 0

    for j in range(1, len(ls)):
    	if ls[j] > (runSum + ls[j]):
            runSum = ls[j]
            i = j
        else:
            runSum += ls[j]

        if runSum > maxSum:
            maxSum = runSum
            start = i
            finish = j

    return start, finish


def pick_peak(above_median_starts, above_median_ends, sub_wps):
    '''
        from region that has 50 < size < 150,
        pick best peak (with maximum wps score)
    '''
    sub_wps = np.asarray(sub_wps)
    above_median_ends = np.asarray(above_median_ends)
    above_median_starts = np.asarray(above_median_starts)
    max_wps_array = np.array([sub_wps[s:e].max() for s, e in izip(above_median_starts, above_median_ends)])
    maximum_wps = np.where(max_wps_array == max_wps_array.max())[0]
    return above_median_starts[maximum_wps], above_median_ends[maximum_wps]


def calling_peaks(chromosome, wpsArray, peak_start, peak_end, peak_count, outFile, peak_size_filter):
    '''
        using peak start and end from wps array,
        find maximum sub array
        and export start and end from maximum subarray
        peak score determine from maximum wps score. 
    '''
    sub_wps = wpsArray[peak_start:peak_end]
    median_sub_wps = np.median(sub_wps)
    adjusted_sub_wps = sub_wps - median_sub_wps
    above_median_starts, above_median_ends = find_peak_region(adjusted_sub_wps)

    if len(above_median_starts)>len(above_median_ends):
        above_median_ends = np.append(above_median_ends,len(adjusted_sub_wps)) 
    if not peak_size_filter:
        above_median_starts, above_median_ends =  pick_peak(above_median_starts, above_median_ends, sub_wps) 

    for above_median_start, above_median_end in izip(above_median_starts, above_median_ends):
        sub_peak_wps = sub_wps[above_median_start:above_median_end]
        nucleosome_start , nucleosome_end = maxSubArray(sub_peak_wps)


        #adjust coordinate
        nucleosome_start, nucleosome_end = peak_start + above_median_start + np.array([nucleosome_start, nucleosome_end])
        nucleosome_center = int((nucleosome_start + nucleosome_end) /2)
        peak_center = (nucleosome_start + nucleosome_end)/2
        nucleosome_size = abs(nucleosome_end - nucleosome_start)
        if (peak_size_filter and 50 < nucleosome_size  < 150 ) or (not peak_size_filter and nucleosome_size > 5):
            peak_score = wpsArray[nucleosome_start:nucleosome_end].max()
            peak_count += 1
            peak_name = '%s_peak%i' %(chromosome, peak_count)
            line = '\t'.join(map(str,[chromosome, nucleosome_start, nucleosome_end, peak_name, peak_score, '+', peak_center]))
            outFile.write(line+'\n')
    return peak_count


def findIntercepts(wpsArray, outFile, chromosome,samplename):
    '''
    looking for the positions that across 0 in the de-noised signal
    '''
    wpsArray = np.asarray(wpsArray)
    start, end = find_peak_region(wpsArray)
    start, end = merge_peaks(start, end)
    peak_count = 0
    #lowPeakBound, upperPeakBound = (15, 120) if 'Long' in lenType else (50,150)
    for peak_start, peak_end  in izip(start, end):
        peak_size = np.abs(peak_end - peak_start)
        if 50 <= peak_size <= 150:
            peak_count = calling_peaks(chromosome, wpsArray, peak_start, peak_end, peak_count, outFile, False)
        elif 150 < peak_size <= 450:
            peak_count = calling_peaks(chromosome, wpsArray, peak_start, peak_end, peak_count, outFile, True)

    printMessage('Written %i peaks to %s' %(peak_count, outFile.name), samplename)
    return 0


def write_peaks(outputWig, outputBed, samplename, lenType):
    bw = pbw.open(outputWig)
    chrom, length = bw.chroms().items()[0]
    chromArray = np.array(bw.values(chrom,0,length))

    with open(outputBed,'w') as outBed:
        if 'Long' in lenType:
            findIntercepts(chromArray, outBed, chrom, samplename)
#        elif 'Short' in lenType:
#
    return 0


def writeWig(chromArray, outputWig, chromosome, samplename):
    outWig =  pbw.open(outputWig,'w')
    chrom_length = len(chromArray)
    outWig.addHeader([(chromosome,chrom_length)])
    outWig.addEntries(chromosome, range(chrom_length), values=map(np.float64,chromArray), span=1)
    outWig.close()
    return 0

def make_wps_array(tempBam, chromosome, tssWindow, wpsWindow, 
                upperBound, lowerBound, lenType, samplename):
    wpsWindow = wpsWindow + 1
    halfWPSwindow = np.divide(wpsWindow,2)
    halfTSSwindow = (tssWindow-1)/2
    with pysam.Samfile(tempBam,'rb') as bam:
        chroms = np.array(bam.references)
        if chromosome not in chroms:
            sys.exit('Wrong chromosome name: %s!' %chromosome)
        chromLength = np.array(bam.lengths)
        chromSize = int(chromLength[chroms==chromosome][0])
        chromArray = extractTSSaln(bam, tssWindow, wpsWindow, halfWPSwindow, upperBound,
                lowerBound, chromosome, chromSize, samplename)
    return chromArray

def runFile(tempBam, outprefix, genome, tssWindow, samplename, chromosome, upperBound, lowerBound, lenType, wpsWindow):
    ''' using bam_file from bed->bam, and analyze the given length type of fragment
        1. Short (35-80bp) WPS_window = 16
        2. Long (120-180bp) WPS_window=120
    '''
    output_folder = os.path.dirname(outprefix)
    bed_folder = output_folder + '/bed_files'
    bigwig_folder = output_folder + '/bigWig_files'
    outputBed = bed_folder + '/'  + samplename + '.'+lenType.split(' ')[0] +'.bed'
    outputWig = bigwig_folder + '/' + samplename +  '.' +lenType.split(' ')[0] +'.bigWig'
    map(make_folder, [bigwig_folder, bed_folder])
    chromArray = make_wps_array(tempBam, chromosome, tssWindow, wpsWindow, 
                                upperBound, lowerBound, lenType, samplename)
    printMessage('Finished calculating WPS for chromosome %s' %(chromosome), samplename)
    writeWig(chromArray, outputWig, chromosome, samplename)
    printMessage('Witten %s' %outputWig, samplename)
    write_peaks(outputWig, outputBed, samplename, lenType)
    printMessage('Witten %s' %outputBed, samplename)
    return 0

def main(inFile, outprefix, genome, tssWindow, chromosome):
    '''
    main function for controling the work flow
    '''
    samplename = os.path.basename(inFile).replace('.bed','')
    printMessage( 'Saving all result to: %s' %outprefix, samplename)
    tempBam = makeBam(inFile, outprefix, genome, samplename)
    run_file = partial(runFile, tempBam, outprefix, genome, tssWindow, samplename  , chromosome)
    upperBound, lowerBound = [80, 180], [35, 120]
    lenType, wpsWindow = ['Short (35-80bp)','Long (120-180bp)'],[16,120]
    map(run_file, upperBound, lowerBound, lenType, wpsWindow)
    return 0

if __name__ == '__main__':
    inFile, outprefix, genome, tssWindow, chromosome = getOpt()
    main(inFile, outprefix, genome, tssWindow, chromosome)
