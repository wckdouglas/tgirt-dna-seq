#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import glob
import os
import re
from multiprocessing import Pool
from scipy.ndimage.interpolation import shift
sns.set_style('white')

def readFile(bedFile):
    basename =  os.path.basename(bedFile)
    samplename = basename.split('.')[0]
    chromosome = basename.split('.')[-3]
    print 'Reading %s for chr: %s ' %(samplename, chromosome)
    column_names = ['chrom','start','end','name','center','strand','intensity']
    df = pd.read_table(bedFile, names=column_names,) \
        .assign(peak_width = lambda d: np.abs(d['end'] - d['start']))\
        .pipe(lambda d: d[(d['peak_width'] >= 50) & (d['peak_width'] <= 500)])\
        .assign(peak_next = lambda d:  np.abs(shift(d['center'],-1) - d['center']))\
        .assign(peak_previous = lambda d: np.abs(shift(d['center'],1) - d['center']))\
        .assign(distance = lambda d: d[['peak_next','peak_previous']].min(axis=1))\
        .assign(samplename = samplename) \
        .pipe(lambda d: d[['distance','samplename']]) \
        .assign(nucleosome_count = 1)\
        .groupby(['distance','samplename']) \
        .agg({'nucleosome_count':np.sum})\
        .reset_index() 
    return df


def normalize_count(d):
    d['normalized_count'] = d.nucleosome_count / d.nucleosome_count.sum()
    return d


def plotDistance(df, figurename):
    upper_bound = 450
    df = df[(df['distance'] <  upper_bound)]  \
        .groupby(['distance','samplename']) \
        .apply(normalized_count)
    with sns.plotting_context('paper',font_scale=1.3):
        p = sns.FacetGrid(data = df, hue='samplename',aspect = 1.6)
        p.map(plt.plot,'distance', 'normalized_count')
        p.set(xlim=(50,500),
             xticks=np.arange(50,500,50))
        p.set_xticklabels(rotation=45)
        p.set_xlabels('Distance to the nearest nucleosome call')
        p.set_ylabels('Fraction of nucleosomes')
        p.add_legend()
        p.savefig(figurename)
    print 'Saved %s' %figurename

def main():
    projectPath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    datapath = projectPath + '/genomeWPS/bed_files'
    resultpath = projectPath + '/figures'
    figurename = resultpath + '/peakDistance.pdf'
    tablename = resultpath + '/peakDistance.tsv'
    if not os.path.isdir(resultpath):
        os.mkdir(resultpath)
    bedFiles = np.array(glob.glob(datapath + '/*.Long.bed'),dtype='string')
    chromosomes = np.array(map(lambda x: x.split('.')[-3], bedFiles),dtype='string')
    usable_chromosomes = np.append(['X','Y'],np.array(np.arange(23),dtype='string'))
    bedFiles = bedFiles[(np.in1d(chromosomes,usable_chromosomes))]
    p = Pool(12)
    dfs = map(readFile, bedFiles)
    p.close()
    p.join()
    df = pd.concat(dfs)
    df.to_csv(tablename,sep='\t',index=False)
    df = pd.read_csv(tablename,sep='\t')
    plotDistance(df, figurename)

if __name__ == '__main__':
    main()
