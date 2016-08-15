#!/usr/bin/env python -u

import matplotlib
matplotlib.use('Agg', warn=False)
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
afont = {'fontname':'Arial'}

def getIsize(bamFile, samplename):
    print 'Analyzing %s ' %samplename
    with pysam.Samfile(bamFile, 'rb') as bam:
        isize = np.array([abs(aln.isize) for aln in bam if (abs(aln.isize) < 80) and \
                                                        aln.is_paired and \
                                                        aln.is_read1 and \
                                                        aln.is_proper_pair and \
                                                        not aln.is_secondary and \
                                                        not aln.is_supplementary],dtype=int)
    print 'Extacted %i alignments from %s' %(len(isize), samplename)
    return isize

def plotFigure(figurename, df):
    plt.figure()
    print 'Normalized counts'
    with sns.plotting_context('paper', font_scale=1.7):
        p = sns.FacetGrid(data=df, col = 'enzyme', row = 'substrate', margin_titles=True)
    p.map(plt.plot, 'size', 'normCount')
    p.map(plt.fill_between, 'size', 'normCount')
    [plt.setp(ax.texts, text="") for ax in p.axes.flat]
    p.set_titles(row_template='{row_name}', col_template="{col_name}", 
            fontweight='bold', size=20, **afont)
    p.set_xticklabels(rotation=90)
    p.fig.text(x = 0.4, y = 0, s = 'Insert Size', fontsize=17, **afont)
    p.fig.text(x = 0, y = 0.7, s = 'Percent total reads', rotation=90, fontsize=17, **afont)
    p.set_xlabels(' ')
    p.set_ylabels(' ')
    p.savefig(figurename)
    print 'Saved %s' %figurename

def bamToDF(bamFile):
    samplename = os.path.basename(bamFile).split('.')[0].split('_')[0]
    enzyme = 'TeI4c' if samplename[0] == 'T'  else 'GsI-IIc'
    substrate = 'DNA' if samplename[1] == 'D'  else 'RNA'
    isize = getIsize(bamFile, samplename)
    size, count = np.unique(isize, return_counts=True)
    df = pd.DataFrame({'size':size,'count':count})
    df['samplename'] = np.repeat(samplename,len(df))
    df['enzyme'] = np.repeat(enzyme, len(df))
    df['substrate'] = np.repeat(substrate, len(df))
    print 'Finished %s ' %samplename
    return df

def normCount(df):
    df['normCount'] = df['count'] / np.sum(df['count']) * 100
    return df

def main():
    start = time.time()
    bampath = '/stor/work/Lambowitz/cdw2854/processivity/filtered_bams'
    figurepath = bampath
    bamFiles = glob.glob(bampath + '/*.bam')
    if len(bamFiles) == 0:
        sys.exit('No bam files found!')
    dfs = Pool(24).map(bamToDF, bamFiles)
    df = pd.concat(dfs)\
            .groupby(['enzyme','substrate','size'])\
            .agg({'count':'sum'})\
            .reset_index() \
            .groupby(['enzyme','substrate']) \
            .apply(normCount)\
            .reset_index()
    tablename = bampath + '/isizeTable.tsv'
    df.to_csv(tablename, sep='\t', index=False)
    figurename = figurepath + 'isizeTable.png'
    plotFigure(figurename, df)
    print 'Time lapsed: %.3f sec' %(time.time() - start)

if __name__ == '__main__':
    main()
