#!/usr/bin/env python -u

from matplotlib import use as mpl_use
mpl_use('Agg')
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
from functools import partial
from collections import defaultdict
sns.set_style('white')

def get_isize(bamFile, samplename):
    print 'Analyzing %s ' %samplename
    isize_dict = defaultdict(lambda: defaultdict(int))
    aln_count = 0
    with pysam.Samfile(bamFile, 'rb') as bam:
        for aln in bam:
            if (abs(aln.isize) <= 500 and abs(aln.isize) >= 10):
                isize_dict[aln.reference_name][abs(aln.isize)] += 1
                aln_count += 1
    print 'Extacted %i alignments from %s' %(aln_count, samplename)
    return isize_dict

def plot_figure(figurename, df):
    with sns.plotting_context('paper',font_scale=1.5):
        p = sns.FacetGrid(data=df.sort_values('isize'),
                          col = 'chrom_annot',
                          col_wrap = 2,
                          sharey=False,
                          legend_out = True,
                          hue = 'samplename')
    p.map(plt.plot, 'isize', 'normalized_count',
          alpha=0.5)
    p.add_legend()
    p.set_titles('{col_name}')
    p.savefig(figurename)
    print 'Saved %s' %figurename
    return 0


def norm_count(d):
    d['normalized_count'] = d['count'] / d['count'].sum()
    return d


def assign_chrom(chrom):
    return 'chrM' if chrom == 'MT' else 'Autosomal and Sex Chromosome'


def make_isize_df(isize_dict):
    dfs = []
    for chrom, chrom_isize in isize_dict.iteritems():
        dfs.append( pd.DataFrame({'isize':chrom_isize.keys(),
                      'count':chrom_isize.values()}) \
                   .assign(chrom = chrom))
    df = pd.concat(dfs, axis=0)
    return df


def parse_bam(table_path, bamFile):
    samplename = os.path.basename(bamFile).split('.')[0]
    insert_table_name = table_path + '/' + samplename + '.csv'
    isize_dict = get_isize(bamFile, samplename)
    df = make_isize_df(isize_dict) \
        .assign(samplename = samplename)
    df.to_csv(insert_table_name, index=False)
    print 'Finished %s and written %s ' %(samplename,insert_table_name)
    return insert_table_name

def makedir(directory):
    os.system('mkdir -p %s' %directory)

def make_joint_table(df_names):
    df = pd.concat(map(pd.read_csv, df_names),
                   axis=0) \
        .reset_index() \
        .drop(['index'],axis=1)\
        .assign(chrom_annot = lambda d: map(assign_chrom, d.chrom)) \
        .groupby(['chrom_annot','samplename']) \
        .apply(norm_count)
    print 'Normalized counts'
    return df


def main():
    start = time.time()
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    data_path = project_path + '/bamFiles'
    table_path = project_path + '/insertSize'
    figure_path = project_path + '/figures'
    figurename = figure_path + '/isizeTable.png'
    tablename = table_path + '/isizeTable.tsv'
    map(makedir, [figure_path, table_path])
    bamFiles = glob.glob(data_path + '/*.bam')
    pool = Pool(24)
    parseBam = partial(parse_bam, table_path)
    df_names = pool.map(parseBam, bamFiles)
    pool.close()
    pool.join()
    df = make_joint_table(df_names)
    df.to_csv(tablename, sep='\t', index=False)
    print 'Written %s' %tablename
    df = pd.read_table(tablename,sep='\t')
    plot_figure(figurename, df)
    print 'Time lapsed: %.3f sec' %(time.time() - start)

if __name__ == '__main__':
    main()
