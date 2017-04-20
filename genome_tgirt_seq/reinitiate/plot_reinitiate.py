#!/usr/bin/env python

from matplotlib import use as mpl_use, ticker
mpl_use('Agg')
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import poisson, chisquare
import glob
import os
import seaborn as sns
import numpy as np
sns.set_style('white')

padding = pd.DataFrame({'fragment_counts':range(5,101), 'counts':[0]*96})
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
    return 'Simulation (reads)' if 'sim' in x else 'TGIRT-seq (UMI)'

def make_simulation_data():
    sim = pd.DataFrame({
            'counts':[18511210, 132626, 2130, 51, 0, 0],
            'fragment_counts':[1 ,2, 3,4, 5, 6]
            }) \
        .append(padding) \
        .assign(normalized_count = lambda s: s['counts']/np.sum(s['counts'])) \
        .assign(samplename ='simulation')
    return sim


def make_cdf(df):
    df['cdf_count'] = np.cumsum(df['normalized_count'])
    return df


def plot_bar(d, k, pv, figurename):
    plt.figure()
    fs=20
    palette = sns.color_palette(['black','green'])
    with sns.plotting_context('paper', font_scale=2):
        ax = sns.barplot(data = d,
                     x = 'fragment_counts',
                     y = 'normalized_count',
                     hue='samplename',
                     palette=palette)
    ax.annotate('$\chi^{2}$: %.3f\nP-value: %.3f' %(k,pv),
                xy=(2,1), fontsize=fs)
    ax.legend(title= ' ', fontsize=fs, loc = (0.5,0.4))
    ax.set_xlim(-0.5,5)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.set_yscale('log')
    ax.set_xlabel('Count per unique fragment',fontsize=fs)
    ax.set_ylabel('% Fragments', fontsize=fs)
    #x = ax.set_xticklabels(ax.get_xmajorticklabels(),rotation=0)
    ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda y,pos: ('{{:.{:1d}f}}'\
                                                                     .format(int(np.maximum(-np.log10(y),0))))\
                                                      .format(y)))

    #ax.tick_params(axis='both', which='major', labelsize=fs)
    #ax.ticklabel_format(style='plain')
    plt.savefig(figurename, transparent=True,bbox_inches='tight')
    print 'Plotted %s' %figurename
    return 0


def main():
    data_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/umi2id_30/frag_counts'
    figurename = data_path + '/reinitiate_model.pdf'
    bar_figurename = figurename.replace('.pdf','_bar.pdf')
    sim_Df = make_simulation_data()
    df = map(read_file, glob.glob(data_path + '/*UMI_1*tsv'))
    df = pd.concat(df, axis = 0)
    df = pd.concat([df,sim_Df],axis=0) \
            .assign(normalized_count = lambda d: d.normalized_count * 100)

    how_many = 5
    chi_df = df.query('fragment_counts < %i' %how_many)
    k, pv = chisquare(chi_df[~chi_df.samplename.str.contains('sim')].normalized_count, 
                      chi_df[chi_df.samplename.str.contains('sim')].normalized_count)

    d = df[df.samplename.str.contains('UMI_1|sim')]\
        .assign(samplename = lambda d: map(rename, d.samplename))

    plot_bar(d, k, pv, bar_figurename)
    return 0

if __name__ == '__main__':
    main()
