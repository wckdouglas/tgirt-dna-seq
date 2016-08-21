#!/bin/env python

# This is the tools for extraciting wps score from bed file
# bam file and calculate window prediction scores from the selected window size
# input: bed file for alignment, each record is a fragment
#        gtf file
#        fai file

import numpy as np
from pybedtools import BedTool
import pysam
import os

def makeBamFilterShortRNA(inFile,outprefix, genome):
    '''
    bed file as input first sort then bedtobam
    index the bam file finally
    '''
    bamFolder = '/'.join(outprefix.split('/')[:-1]) +'/bamFiles'
    tempBam = bamFolder + '/' + outprefix.split('/')[-1] + '.filtered.sorted.bam'
    index_name = tempBam + '.bai'
    small_RNA = '/Users/wckdouglas/plasmaDNA/reference/smallRNA.bed'
    small_RNA_bed = Tool(small_RNA)
    if not os.path.isfile(tempBam):
        print 'Making %s ' %tempBam
        BedTool(inFile)\
            .sort()\
            .to_bam(g=genome)\
            .intersect(b=small_RNA_bed, v=True, f=0.8,r=True,s=True) \
            .saveas(tempBam)
        if os.path.isfile(index_name):
            os.remove(index_name)
        index = pysam.index(tempBam)
    else:
        print 'Used existing bamfile: %s' %tempBam
    return tempBam

def writeBam(inFile, tempBam, genome):
    print 'Making %s ' %tempBam
    index_name = tempBam + '.bai'

    temp_dir = tempBam + '.temp'
    os.mkdir(temp_dir)
    command = 'cat %s' %inFile +\
        '| sort -k 1,1 -k2,2n --temporary-directory=%s' %(temp_dir) +\
        '| bedtools bedtobam -i /dev/stdin -g %s' %(genome) +\
        '> %s' %(tempBam)
#    BedTool(inFile)\
#        .sort()\
#        .to_bam(g=genome)\
#        .saveas(tempBam)
    os.system(command)
    os.rmdir(temp_dir)
    print 'Converted to: %s' %tempBam
    if os.path.isfile(index_name):
        os.remove(index_name)
    index = pysam.index(tempBam)
    print 'Indexed: %s' %tempBam
    return 0

def makeBam(inFile, samplename, genome):
    '''
    bed file as input first sort then bedtobam
    index the bam file finally
    '''
    bamFolder = os.path.dirname(inFile) + '/bamFiles'
    tempBam = bamFolder + '/' + samplename + '.sorted.bam'
    if not os.path.isfile(tempBam):
        writeBam(inFile, tempBam, genome)
    else:
        inFileTime = os.path.getmtime(inFile)
        bamFileTime = os.path.getmtime(tempBam)
        if inFileTime > bamFileTime:
            print 'Bed file newer than bam %s' %tempBam
            writeBam(inFile, tempBam, genome)
        else:
            print 'Used existing bamfile: %s' %tempBam
    return tempBam

def pushWPStoArray(aln, halfWPSwindow, posStart, posEnd, tssWindow, tss, isize, wpsWindow):
    """
    for a given alignment, compute the regions that can be fully aligned and not
    e.g. [-1, -1 , -1, 1 , 1, 1, 1, 1, 1, -1, -1, -1] for a wps window -f 6 (halfwindow 3 )
    this will be added to the defined transcription start site wps score array after
    adjusting fot the position
    """
    transcriptAlnWPS = np.zeros(tssWindow) # setting the tss window as zeros wps array
    alnWPS = np.zeros(isize + wpsWindow) #adding halfwindow to both side of the alignment
    alnWPS[wpsWindow:-wpsWindow] = 1 # only half window after aln start, the alignment is fully covering the wps window
                                     # and we added half window on the previous line
    alnWPS[alnWPS != 1 ] = -1 #making the alignment wps with ends containing windows
    alnShift = posStart - (aln.pos - halfWPSwindow) # the distance between alignment start and right side of the wps window:
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

def calculateWPS(bam, tss, tssWindow, wpsWindow, halfWPSwindow, upperBound, lowerBound):
    '''
    for each gene start site region:
        for i in each position:
            fetch all alignments in the region with the window size (e.g. 120)
            calculate number of alignments fully span and partially mapped
            get wps[i] = number of full map - number of end mapped
    '''
    transcriptWps = np.zeros(tssWindow)
    posStart = int(tss.start)  - 1 # accommodate for pysam 0-base coordinate
    posEnd = int(tss.end)  - 1
    coverage = 0
    for aln in bam.fetch(reference = str(tss.chrom), start = (posStart - halfWPSwindow), end = (posEnd+halfWPSwindow)):
        isize = abs(aln.qlen)
        if upperBound > isize > lowerBound:
            transcriptWps += pushWPStoArray(aln, halfWPSwindow, posStart, posEnd, tssWindow, tss, isize, wpsWindow)
            coverage += 1
    transcriptWps = np.flipud(transcriptWps) if str(tss.strand) == '-' else transcriptWps
    return transcriptWps, coverage
