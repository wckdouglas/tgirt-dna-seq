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

data_path = '/stor/work/Lambowitz/cdw2854/genomeDNA/clustered_map/rmdup_bam'
df = pd.read_table(data_path + '/bad_promoter_coverage.tsv')
df = df \
    .assign(seq_type = lambda d: map(assign_seq_type,d['samplename'])) \
    .assign(position = lambda d: d['position']-1900) \
    .groupby(['seq_type', 'position']) \
    .agg({'coverage':np.sum}) \
    .reset_index() \
    .groupby('seq_type') \
    .apply(norm_df)

sns.set_style('white')
with sns.plotting_context('paper',font_scale = 2):
    p = sns.FacetGrid(data = df, hue = 'seq_type', 
                    legend_out = True, aspect = 2) 
p.map(plt.plot, 'position', 'coverage')
sns.plt.axvline(x=0, ymin = 0, ymax = 10000, c='red')
sns.plt.axvline(x=200, ymin = 0, ymax = 10000, c='red')
sns.plt.suptitle('1000 Bad promoter region\n(Aggregated coverage)', 
                fontweight='bold', fontsize=18,y=1.08)
p.add_legend(title=' ')
p.set_xlabels('Position')
p.set_ylabels('Normalized coverage') 
p.set_xticklabels(rotation = 75)
figurename = data_path + '/bad_promoter_comparison.pdf'
p.savefig(figurename)
print 'Plotted ', figurename
