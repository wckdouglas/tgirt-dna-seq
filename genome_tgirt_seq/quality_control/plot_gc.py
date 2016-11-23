#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import pandas as pd
import os
import numpy as np
import glob
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('white')

def rename_function(x):
    if 'SRR' in x:
        return 'NEBNext'
    elif 'no_bias' in x: 
        return 'Simulation: No bias'
    elif 'sim_ligation' in x:
        return "Simulation: Only 5' bias"
    elif x == 'EG_sim':
        return 'Simulation: Both end bias'
    elif 'template' in x:
        return "Simulation: Only 3' bias"
    elif 'EG1' in x:
        return 'TGIRT-seq'
    elif 'EG2' in x:
        return 'TGIRT-seq 2'
    else:
        return x

def normalize(x):
    x = np.array(x)
    return np.true_divide(x, x.sum())

def read_gc_table(filename):
    samplename = os.path.basename(filename).split('.')[0]
    df = pd.read_table(filename, skiprows = 6) \
        .assign(samplename = samplename) \
        .assign(normalize_windows = lambda d: normalize(d['WINDOWS']))
    return df


project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome/'
picard_path = project_path + '/picard_results'
figure_path  = project_path + '/figures'
figurename = figure_path + '/gc_plot.pdf'
table_names = glob.glob(picard_path + '/*.txt')
dfs = map(read_gc_table, table_names)
df = pd.concat(dfs, axis = 0) \
        .assign(samplename = lambda d: map(rename_function, d.samplename))
df.to_csv(figurename.replace('.pdf','.csv'),index=False)
df = df.query('samplename != "TGIRT-seq 2"')

windows = dfs[0]['normalize_windows'].rolling(window=10,center=False).mean() * 10
with sns.plotting_context('paper', font_scale = 2):
    p = sns.FacetGrid(data = df, hue = 'samplename', 
            legend_out=False, aspect = 1.3, size = 6)
p.map(plt.plot,'GC','NORMALIZED_COVERAGE', alpha=0.7)
plt.plot( np.arange(len(windows)), windows, color = 'salmon')
plt.fill_between(np.arange(len(windows)),y1=0, y2 = windows, color = 'salmon', alpha = 0.5)
plt.axhline(y = 1, xmin = 0, xmax = 100, linestyle='--', alpha = 0.7, color = 'grey')
p.set(xlabel = 'GC %', ylabel = 'Normalized Coverage')
p.add_legend(title = ' ')
p.savefig(figurename)
print 'Saved: %s' %figurename
