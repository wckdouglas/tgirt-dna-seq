#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
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
    columnNames = ['chrom','start','end','name','center','strand']
    df = pd.read_csv(bedFile,sep='\t', names=columnNames)
    df['nucleosome'] = np.abs(df['end'] - df['start'])
    df = df[(df['nucleosome'] >= 50) & (df['nucleosome'] <= 500)]
    df['peakNext'] = np.abs(shift(df.center,-1) - df.center)
    df['peakPrevious'] = np.abs(shift(df.center,1) - df.center)
    df['distance'] = df[['peakNext','peakPrevious']].min(axis=1)
    matrix = np.array([df['distance'], np.repeat(samplename,len(df))])
    df = pd.DataFrame(matrix.transpose(), columns = ['distance','samplename'])
    return df

def plotDistance(df, figurename):
    lowerBound = 450
    df['distance'] = np.asarray(df['distance'], dtype=int)
    df = df[(df['distance'] <  lowerBound)] 
    df.columns = ['distance','Sample']
    df['plotTitle'] = np.repeat('Nearest nucleosome call',len(df))
    with sns.plotting_context('paper',font_scale=1.3):
        p = sns.FacetGrid(data = df, hue='Sample',aspect = 1.6, legend_out=False)
        p.map(sns.distplot,'distance',hist=False,bins=1)
        p.set(xlim=(50,500),
             xticks=np.arange(50,500,50))
        p.set_xticklabels(rotation=45)
        p.set_xlabels('Distance to the nearest nucleosome call')
        p.set_ylabels('Fraction of nucleosomes')
        p.add_legend()
        p.savefig(figurename)
    print 'Saved %s' %figurename

def main():
    projectPath = '/scratch/cdw2854/plasmaDNA'
    datapath = projectPath + '/genomeWPS'
    resultpath = projectPath + '/figures'
    figurename = resultpath + '/peakDistance.pdf'
    tablename = resultpath + '/peakDistance.tsv'
    if not os.path.isdir(resultpath):
        os.mkdir(resultpath)
    bedFiles = np.array(glob.glob(datapath + '/*.Long.bed'),dtype='string')
    chromosomes = np.array(map(lambda x: x.split('.')[-3], bedFiles),dtype='string')
    bedFiles = bedFiles[(np.in1d(chromosomes,np.arange(23)))]
    bedFiles = filter(lambda x: re.search('SRR|NT|PD|RNa',x), bedFiles )
    dfs = Pool(24).map_async(readFile, bedFiles).get()
    df = pd.concat(dfs)
    df.to_csv(tablename,sep='\t',index=False)
    df = pd.read_csv(tablename,sep='\t')
    plotDistance(df, figurename)

if __name__ == '__main__':
    main()
