#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from math import sqrt
sns.set_style('white')

number_of_promoter = 100
tablename = 'bad_promoter.tsv'
figurename = 'promoter.pdf'
df = pd.read_table('bad_promoter_coverage.tsv')

aggregated = df \
    .query("name != 'TMEM163'")\
    .groupby(['name']) \
    .agg({'coverage':np.sum}) \
    .reset_index()\
    .sort_values(['coverage'],ascending = False)\
    .head(number_of_promoter)\
    .drop(['coverage'],axis=1)\
    .merge(df, on = 'name', how = 'inner') \
    .assign(name = lambda d: map(lambda x, y: x + ': '+str(y), d['name'], d['ratio']))

sorted_name = aggregated \
    .pipe(lambda d: d[['ratio','name']]) \
    .drop_duplicates()\
    .sort_values(['ratio'], ascending=True)
    
with sns.plotting_context('paper'):
    p = sns.FacetGrid(data = aggregated, col = 'name',
        sharey=False,
        col_wrap = int(sqrt(number_of_promoter)),
        col_order = sorted_name['name'])
p.map(plt.plot,'position','coverage')
p.map(sns.plt.axvline, x = 1900, c ='salmon')
p.map(sns.plt.axvline, x = 2100, c ='salmon')
p.set_titles('{col_name}')
p.savefig(figurename)
print 'Saved:', figurename

figurename = 'aggregated_promoter.pdf'
aggregated = df \
    .pipe(lambda d: d[d['ratio']<0.2])\
    .groupby(['position']) \
    .agg({'coverage':np.sum})  \
    .reset_index()

with sns.plotting_context('paper'):
    p = sns.FacetGrid(data = aggregated)
p.map(plt.plot,'position','coverage')
p.map(sns.plt.axvline, x = 1900, c ='salmon')
p.map(sns.plt.axvline, x = 2100, c ='salmon')
p.savefig(figurename)
print 'Saved:', figurename
