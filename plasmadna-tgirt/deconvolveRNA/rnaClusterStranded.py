#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pybedtools import BedTool
import glob
import os
import pandas as pd
import numpy as np
import pylab as plt
from multiprocessing import Pool
import seaborn as sns
from dmisc import folder_action
import rnaCluster  
sns.set_style('white')

def findRNA(args):
    bedfile, rnaBed, motifBed = args
    samplename = os.path.basename(bedfile).split('.')[0]
    print 'Running %s' %samplename
    shortBed = BedTool(bedfile)
    longBed = BedTool(bedfile.replace('Short','Long'))
    overlap_fraction = 0.000001
    nonDNA = shortBed.intersect(b=motifBed, f=0.8, r = True, v=True) \
            .intersect(b=longBed, v=True)
    senseIntersection= nonDNA.intersect(b=rnaBed, f = overlap_fraction, wb=True, s=True) 
    antisenseIntersection = nonDNA.intersect(b=rnaBed, f = overlap_fraction, wb=True, S=True) 
    names = map(lambda x: samplename + '_' + x, ['sense','antisense'])
    dfs = map(rnaCluster.intersect2dataframe,
            [senseIntersection, antisenseIntersection],
            names)
    df = pd.concat(dfs, axis = 0)
    return df

def main():
    projectpath = '/scratch/02727/cdw2854/plasmaDNA'
    bedpath = projectpath + '/genomeWPS/mergedStrandedBED'
    resultpath = projectpath + '/rnaCallsStrand'
    referencepath = '/corral-repl/utexas/2013lambowitz/Ref/GRCh38/Bed_for_counts_only'
    bedfiles = glob.glob(bedpath + '/*bed')
    rnaRef = referencepath + '/genes.bed'
    motifRef = referencepath + '/MotifFeatures.gff'
    rnaBed = BedTool(rnaRef)
    motifBed = BedTool(motifRef)
    resultFilename = resultpath + '/rnaCount_strand.tsv'
    figurename = resultpath + '/rnaCount_strand.png'
    pool = Pool(24)
    folder_action.makeFolder(resultpath)
    dfs = pool.map(findRNA, [(bedfile, rnaBed, motifBed) for bedfile in bedfiles])
    df = pd.concat(dfs)
    df.to_csv(resultFilename, sep='\t', index=False)
    print 'Written %s' %resultFilename
    rnaCluster.plotType(df, figurename)
    return 0

if __name__ == '__main__':
    main()
