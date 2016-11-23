#!/usr/bin/env python

import matplotlib 
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np 

def assign_seq_type(x):
    return 'Illumina' if 'illumina' in x else 'TGIRT'

def norm_df(d):
    d['coverage'] = np.true_divide(d['coverage'],np.max(d['coverage']))
    return d

data_path = '/Users/wckdouglas/Desktop'
df = pd.read_table(data_path + '/bad_promoter_coverage.tsv') \
    .assign(seq_type = lambda d: map(assign_seq_type,d['samplename'])) \
    .groupby(['name','position','seq_type']) \
    .agg({'coverage':np.sum}) \
    .reset_index()
    
df = df \
    .groupby(['name']) \
    .agg({'coverage':np.mean}) \
    .reset_index() \
    .sort_values(['coverage'], ascending = False) \
    .head(100) \
    .pipe(lambda d: d[['name']]) \
    .merge(df, on = 'name', how = 'inner')\
    .groupby(['seq_type','name']) \
    .apply(norm_df) \
    .assign(position = lambda d: d['position'] - 1900)

p = sns.FacetGrid(data = df, col = 'name', 
                hue = 'seq_type', col_wrap = 10)
p.map(sns.plt.plot, 'position', 'coverage')
p.add_legend()
for ax in p,fig.axes:
    ax.axvline(x= 0, ymin = 0, ymax = 10000, c='red')
    ax.axvline(x= 200, ymin = 0, ymax = 10000, c='red')

