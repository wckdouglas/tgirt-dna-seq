#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import numpy as np
import seaborn as sns
import fileinput
import pylab as plt
from pybedtools import BedTool
import os
import pandas as pd
import time
from multiprocessing import Pool
from functools import partial
sns.set_style('white')

def plottingFigure(figurename, df):
    df.columns = ['distance','Chromosomes'] 
    with sns.plotting_context('paper',font_scale=1.6):
        p = sns.FacetGrid(data = df, hue = 'Chromosomes', aspect = 2, size = 7, legend_out=True)
        p.map(sns.distplot, 'distance', hist=False)
        p.set(xlim=(-1000,1000),
                xlabel = 'Distance from NT to closest SRR2130051 nucleosom calls',
                ylabel ='Fraction of nucleosomes',
                xticks = np.arange(-1100,1100,100))
        p.set_xticklabels(np.arange(-1100,1100,100),rotation=45)
        p.add_legend()
        plt.savefig(figurename)
    print 'Plotted %s' %figurename
    return 0

def parseClosetLine(closestLine):
    fields = closestLine.fields
    center1 = int(fields[4])
    center2 = int(fields[11])
    distance = center2 - center1
    return distance

def closestPeak(bedpath, file1, file2, chromosome):
    print 'Running chromosome: %s' %chromosome
    bedFile1 = '%s/%s.%s.Long.bed' %(bedpath, file1, chromosome)
    bedFile2 = bedFile1.replace(file1,file2)
    if os.path.isfile(bedFile1) and os.path.isfile(bedFile2):
        bed1 = BedTool(bedFile1)
        bed2 = BedTool(bedFile2)
        centerDistance = np.array([parseClosetLine(line) for line in BedTool(bed1).closest(bed2)])
        centerDistance = centerDistance[(-1000 < centerDistance) & (centerDistance< 1000)]
        df = pd.DataFrame(centerDistance, columns=['distance'])
        df['chrom'] = chromosome
        return df

def makeDir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    return 0

def main():
    start = time.time()
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bedpath = projectpath + '/genomeWPS'
    figurepath = projectpath + '/figures'
    figurename = figurepath + '/predictedNucleosomeDistance.pdf'
    tablename = figurename.replace('pdf','tsv')
    makeDir(figurepath)
    file1 = 'PD-merged'
    file2 = 'SRR2130051'
    chromosomes = map(str,np.arange(1,23))
    chromosomes = np.concatenate([chromosomes,['X','Y']])
    closestPeakFunc = partial(closestPeak, bedpath, file1, file2)
    p = Pool(24)
    dfs = p.map(closestPeakFunc, chromosomes)
    dfs = [df for df in dfs if df is not None]
    df = pd.concat(dfs)
    df.to_csv(tablename, index=False)
    plottingFigure(figurename, df)
    print 'Writtten %s in %.3f min' %(tablename, np.true_divide(time.time() - start,60))
    return 0

if __name__ == '__main__':
    main()
