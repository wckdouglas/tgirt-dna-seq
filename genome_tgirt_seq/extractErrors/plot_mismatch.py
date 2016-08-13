#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
sns.set_style('white')


def plotting(gdf,figurename):
    gdf= gdf[gdf['new_ref'] != gdf['read']]
    with sns.plotting_context('paper',font_scale=1.4):
        p = sns.FacetGrid(data = gdf, col='samplename', size = 7)
    p.map(sns.barplot,'mismatch','count',color='skyblue')
    p.set_xticklabels(rotation=60)
    p.set(xlabel=' ', ylabel='Frequency of errors')
    for ax in p.fig.axes:
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
    p.set_titles('cluster member > {col_name}')
    plt.ticklabel_format(axis='y',style='sci',scilimits=(1,2))
    p.savefig(figurename)
    print 'Saved %s ' % figurename


def normCount(df):
    df['count'] = np.true_divide(df['count'],np.sum(df['count']))
    return df

datapath = '/scratch/02727/cdw2854/jurkatCells/mismatches'
df = pd.read_csv(datapath + '/mismatch.tsv',sep='\t')
df.groupby(['samplename']).agg({'count':sum})
df['mismatch'] = map(lambda x,y: y+' to '+x, df['read'],df['ref'])
df = df.groupby(['samplename']).apply(normCount)
df = df[(df['read'] != df['ref']) & ~(df['mismatch'].str.contains('N'))]
