#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import os
import glob
from multiprocessing import Pool
import pandas as pd
from functools import partial
sns.set_style('white')

def count(f):
    line_count = sum(1 for l in open(f, 'r'))
    return line_count, os.path.basename(f).split('.')[0].replace('_unique','')

def files_to_df(demultiplexed_path, file_type):
    files = glob.glob(demultiplexed_path + '/*%s*.bed' %(file_type))
    count_record = Pool(12).map(count, files) 
    df = pd.DataFrame(count_record, columns = ['fragment_count','name']) \
            .groupby(['name']) \
            .sum() \
            .reset_index() \
            .assign(file_type = file_type)
    return df

def fix_dedup(x):
    return 'Fragment Ends' if x == 'collapse' else 'UMI + Fragment Ends'

def plotting(df, figurename):
    with sns.plotting_context('paper',font_scale = 2.5):
        ax = sns.barplot(data = df, x = 'name', y = 'fragment_count', hue = 'Deduplication')
        plt.legend(title = ' ', loc = (0.4,0.8))
        plt.ylabel('Fragment Count')
        plt.xlabel(' ')
        ax.set_xticklabels(map(lambda x: 'Sample %i' %x, range(1,4)))
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        plt.savefig(figurename, transparent=True)

def main():
    demultiplexed_path = os.environ['WORK'] + '/cdw2854/plasmaDNA/bedFiles/demultiplexed'
    figurename = demultiplexed_path + '/duplicate_count.pdf'
    file_types = ['unique','collapse']
    ftd_func = partial(files_to_df, demultiplexed_path)
    df = pd.concat(map(ftd_func, file_types), axis=0) \
        .pipe(lambda d: d[~d.name.str.contains('1203')]) \
        .assign(Deduplication = lambda d: map(fix_dedup, d['file_type']))
    plotting(df, figurename)
    print 'Plotted %s' %figurename

if __name__ == '__main__':
    main()
