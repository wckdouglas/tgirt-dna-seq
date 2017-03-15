#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson, ks_2samp
import glob
import os
import seaborn as sns
import numpy as np

padding = pd.DataFrame({'fragment_counts':range(4,101), 'counts':[0]*97})
def read_file(filename):
    samplename = os.path.basename(filename.split('.')[0])
    df = pd.read_table(filename, sep=' ',
                       names = ['fragment_counts','counts'])\
        .append(padding) \
        .groupby(['fragment_counts'])\
        .max()\
        .assign(samplename = samplename)\
        .assign(normalized_count = lambda d: d['counts']/np.sum(d['counts']))\
        .reset_index()
    return df

def rename(x):
    return 'Simulation' if 'sim' in x else 'TGIRT-seq'

def make_simulation_data():
    sim = pd.DataFrame({
            'counts':[1264781, 3120, 12],
            'fragment_counts':[1 ,2, 3]
            }) \
        .append(padding) \
        .assign(normalized_count = lambda s: s['counts']/np.sum(s['counts'])) \
        .assign(samplename ='simulation')
    return sim


def make_cdf(df):
    df['cdf_count'] = np.cumsum(df['normalized_count'])
    return df


def plot_line(d, k, pv, figurename):
    d = d.reset_index()\
        .groupby(['samplename']) \
        .apply(make_cdf)
    plt.figure()
    sns.set_style('white')
    with sns.plotting_context('paper',font_scale=1.2):
        p = sns.FacetGrid(data = d, legend_out = False,
                  hue ='samplename')
    p.map(plt.plot, 'fragment_counts','cdf_count', alpha=0.5)
    p.add_legend(title = ' ')
    plt.xlabel('cDNA counts per substrate', fontweight='bold')
    plt.ylabel('Cumulative probability', fontweight='bold')
    p.set(xticks = range(0,6), xlim=(0,5))
    p.fig.axes[0].annotate('KS-stat: %.3f\nP-value: %.3f' %(k,pv), xy=(2,0.996))
    #p.fig.axes[0].spines['top'].set_visible(False)
    #p.fig.axes[0].spines['right'].set_visible(False)
    plt.savefig(figurename, transparent=True)
    print 'Plotted %s' %figurename
    return 0


def plot_bar(d, k, pv, figurename):
    plt.figure()
    fs=20
    ax = sns.barplot(data = d,
                     x = 'fragment_counts',
                     y = 'normalized_count',
                     hue='samplename')
    ax.annotate('KS-stat: %.3f\nP-value: %.3f' %(k,pv),
                xy=(2,0.3), fontsize=fs)
    ax.legend(title= ' ', fontsize=fs, loc = (0.5,0.2))
    ax.set_xlim(-0.5,5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.yscale('log')
    ax.set_xlabel('cDNA counts per substrate',fontsize=fs)
    ax.set_ylabel('Probability', fontsize=fs)
    x = ax.set_xticklabels(ax.get_xmajorticklabels(),rotation=0)
    ax.tick_params(axis='both', which='major', labelsize=fs)
    plt.savefig(figurename, transparent=True,bbox_inches='tight')
    print 'Plotted %s' %figurename
    return 0


def main():
    data_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id_30/frag_counts'
    figurename = data_path + '/reinitiate_model.pdf'
    bar_figurename = figurename.replace('.pdf','_bar.pdf')
    sim_Df = make_simulation_data()
    df = map(read_file, glob.glob(data_path + '/*tsv'))
    df = pd.concat(df, axis = 0)
    df = pd.concat([df,sim_Df],axis=0)

    k, pv = ks_2samp(df[df.samplename.str.contains('kh_1')]['normalized_count'],
               sim_Df.normalized_count)

    d = df[df.samplename.str.contains('kh_1|sim')]\
        .assign(samplename = lambda d: map(rename, d.samplename))

    plot_line(d, k, pv, figurename)
    plot_bar(d, k, pv, bar_figurename)
    return 0

if __name__ == '__main__':
    main()
