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
import re
from scipy.signal import argrelextrema
from multiprocessing import Pool
sns.set_style('white')

def getIsize(bamFile, samplename):
    print 'Analyzing %s ' %samplename
    with pysam.Samfile(bamFile, 'rb') as bam:
        isize_array = [(abs(aln.isize), aln.reference_name) for aln in bam if (abs(aln.isize) < 300 and abs(aln.isize) > 30)]
    print 'Extacted %i alignments from %s' %(len(isize_array), samplename)
    return isize_array

def plotFigure(figurename, df):
    x_axis_range = np.arange(0,301,50)
    with sns.plotting_context('paper',font_scale=1.5):
        p = sns.FacetGrid(data=df, col = 'samplename', col_wrap = 4)
    p.map(plt.plot, 'isize', 'normCount')
    p.map(plt.fill_between, 'isize','normCount')
    p.set_titles('{col_name}')
    p.set_xticklabels(rotation=60)
    p.set(xlabel ='Insert Size', ylabel = 'Relative abundance')
    p.savefig(figurename)
    print 'Saved %s' %figurename
    return 0

def normCount(d):
    d['normCount'] = d['count'] / np.sum(d['count'])
    return d

def renameChrom(x):
    return 'chrM' if x == 'MT' else 'Autosomal and Sex Chromosome'

def parseBam(bamFile):
    samplename = os.path.basename(bamFile).split('.')[0]
    isize_array = getIsize(bamFile, samplename)
    df = pd.DataFrame(isize_array, columns = {'isize','chrom'})\
        .assign(count = 1)\
        .groupby(['isize','chrom'])\
        .agg({'count':np.sum})\
        .reset_index() \
        .assign(chrom = lambda d: map(renameChrom, d['chrom'] )) \
        .assign(samplename = samplename)
    print 'Finished %s ' %samplename
    return df

def makedir(directory):
    os.system('mkdir -p %s' %directory)

def main():
    start = time.time()
    project_path = '/stor/work/Lambowitz/cdw2854/bisufite_seq'
    data_path = project_path + '/rmdupBAM'
    table_path = project_path + '/insertSize'
    figure_path = project_path + '/figures'
    figurename = figure_path + '/isizeTable.png'
    tablename = table_path + '/isizeTable.tsv'
    map(makedir, [figure_path, table_path])
    bamFiles = glob.glob(data_path + '/*.bam')
    #pool = Pool(12)
    #dfs = pool.map(parseBam, bamFiles)
    #pool.close()
    #pool.join()
    #df = pd.concat(dfs)
    #df.to_csv(tablename, sep='\t', index=False)
    df = pd.read_table(tablename,sep='\t')\
        .drop(['chrom'],axis=1) \
        .groupby(['samplename','isize']) \
        .agg({'count':np.sum}) \
        .reset_index() \
        .groupby(['samplename']) \
        .apply(normCount) \
        .reset_index()
    plotFigure(figurename, df)
    print 'Time lapsed: %.3f sec' %(time.time() - start)

if __name__ == '__main__':
    main()
