#!/usr/bin/env python -u

import matplotlib
matplotlib.use('Agg')
import pysam
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import time
import glob
import pandas as pd
import os
from scipy.signal import argrelextrema
from multiprocessing import Pool
sns.set_style('white')

def getIsize(bamFile, samplename):
    print 'Analyzing %s ' %samplename
    with pysam.Samfile(bamFile, 'rb') as bam:
        isize_array = [abs(aln.isize) for aln in bam if (30 < abs(aln.isize) < 300 \
                                                                                and aln.is_proper_pair \
                                                                                and aln.is_read1 \
                                                                                and not aln.is_supplementary \
                                                                                and not aln.is_secondary \
                                                                                and not aln.is_unmapped)]
    print 'Extacted %i alignments from %s' %(len(isize_array), samplename)
    return isize_array

def plotFigure(figurename, df):
    with sns.plotting_context('paper',font_scale=1.5):
        plt.figure()
        print 'Normalized counts'
        p = sns.FacetGrid(data=df, col = 'samplename', col_wrap = 4)
        p.map(plt.plot, 'isize', 'normCount')
        p.set_titles('{col_name}')
        p.savefig(figurename)
        print 'Saved %s' %figurename
    return 0

def normCount(d):
    d['normCount'] = d['count'] / np.max(d['count'])
    return d

def main(bamFile):
    samplename = os.path.basename(bamFile).split('.')[0].split('_')[0]
    print 'Running %s ' %samplename
    isize_array = getIsize(bamFile, samplename)
    df = pd.DataFrame({'isize':isize_array})\
	.assign(count = 1) \
	.groupby(['isize']) \
	.agg({'count':np.sum})\
	.reset_index()  \
	.assign(samplename = samplename) 
    print 'Finished %s ' %samplename
    return df

if __name__ == '__main__':
    start = time.time()
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    datapath = project_path + '/bamFiles'
    figurepath = project_path + '/figures'
    figurename = figurepath + '/isizeTable.png'
    tablename = datapath + '/isizeTable.tsv'
    bamFiles = glob.glob(datapath + '/*bam')
    print datapath
    dfs = Pool(12).map(main, bamFiles)
    df = pd.concat(dfs)
    df.to_csv(tablename, sep='\t', index_label = False, index=False)
    df = pd.read_table(tablename,sep='\t')\
	.groupby(['samplename']) \
	.apply(normCount) \
	.reset_index() 
    plotFigure(figurename, df)
    print 'Time lapsed: %.3f sec' %(time.time() - start)
