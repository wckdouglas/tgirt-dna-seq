#!/usr/bin/env python

from matplotlib import use as mpl_use 
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style('white')

figure_path = '/stor/work/Lambowitz/cdw2854/bisufite_seq/figures'
datapath = '/stor/work/Lambowitz/cdw2854/bisufite_seq/insertSize'
# BS-seq insert size table
meth_df = pd.read_table(datapath+'/P13B_mix_S1.tsv')\
        .sort_values('isize')\
        .assign(data_type = 'Bisulfite-treated cfDNA')\
        .pipe(lambda d: d[['isize','percentage','data_type']])

# normal cfDNA
plasma_df = pd.read_csv('/stor/work/Lambowitz/cdw2854/plasmaDNA/insertSize/P1022_1_S3_umi2id.csv')\
        .groupby(['isize'], as_index=False)\
        .agg({'count':'sum'})\
        .assign(percentage = lambda d: np.true_divide(d['count'],d['count'].sum()))\
        .sort_values('isize')\
        .assign(data_type = 'Plasma cfDNA')\
        .pipe(lambda d: d[['isize','percentage','data_type']]) 


df = pd.concat([plasma_df,meth_df],axis=0)\
        .sort_values('isize')\
        .query('isize > 25') 
# plot isize
with sns.plotting_context('paper',font_scale = 2):
    p = sns.FacetGrid(data = df, 
                  hue = 'data_type',
                  sharey = False,
                  legend_out = False,
                 size = 5,aspect=1.5) 
p.map(plt.plot, 'isize', 'percentage')
p.add_legend(title = ' ',fontsize=15)
p.set(xlabel = 'Fragment size', 
      ylabel = '% fragments')
p.fig.axes[0].set_yticklabels(p.fig.axes[0].get_yticks() *100)
figurename = figure_path + '/insert_size.pdf'
plt.savefig(figurename)
print 'Saved: ', figurename


# read methylation sites
methyl_plasma = '/stor/work/Lambowitz/cdw2854/bisufite_seq/methyl_calling/P13B_mix_S1_CpG.meth.bedGraph'
plasma_methyl_df = pd.read_table(methyl_plasma,
              names = ['chrom','start','end','methyl_dens'],
              skiprows=1) 
plasma_methyl_df.head()

# plot methylation density
with sns.plotting_context('paper', font_scale=1.5):
    p = sns.distplot(plasma_methyl_df['methyl_dens']*100,
                     hist=True, kde=True, bins=50)
plt.yticks(p.get_yticks(), p.get_yticks()*100)
p.set_xlabel('Methylation density')
p.set_ylabel('Percentage')
p.spines['top'].set_visible(False)
p.spines['right'].set_visible(False)
figurename = figure_path + '/methyl_dens.pdf'
plt.savefig(figurename)
print 'Saved: ', figurename
