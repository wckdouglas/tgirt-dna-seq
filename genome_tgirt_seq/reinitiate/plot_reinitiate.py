#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson
import glob 
import os
import seaborn as sns
import numpy as np

def read_file(filename):
    samplename = os.path.basename(filename.split('.')[0])
    df = pd.read_table(filename, sep=' ', names = ['fragment_counts','counts'])\
            .assign(samplename = samplename)\
            .assign(normalized_count = lambda d: d['counts']/np.sum(d['counts']))
    return df

def rename(x):
    if 'sim' in x:
        return 'Simulation'
    else:
        return 'TGIRT-seq (1)' if 'kh_1' in x else 'TGIRT-seq (2)' 

def make_simulation_data():
    sim = pd.DataFrame({
        'counts':[1264781, 3120, 12],
        'fragment_counts':[1 ,2, 3]})
    sim['normalized_count'] = sim['counts']/np.sum(sim['counts'])
    sim['samplename'] = 'simulation'
    return sim

data_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id_30/frag_counts'
df = map(read_file, glob.glob(data_path + '/*tsv'))
df = pd.concat(df, axis = 0)
df = pd.concat([df,make_simulation_data()],axis=0)

d = df[df.samplename.str.contains('kh|sim')]\
        .assign(samplename = lambda d: map(rename, d.samplename))
p = sns.FacetGrid(data = d, legend_out = False,
                  hue ='samplename')
p.map(plt.plot, 'fragment_counts','normalized_count', alpha=0.5)
p.add_legend(title = ' ')
p.set(xticks = range(1,4))
p.set_axis_labels('cDNA counts per substrate','Probability')
figurename = data_path + '/reinitiate_model.pdf'
p.savefig(figurename)
print 'Plotted %s' %figurename
