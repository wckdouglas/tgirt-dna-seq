#!?usr/bin/env python


from matplotlib import rc, use as mpl_use
mpl_use('Agg')
import pandas as pd
import numpy as np
import seaborn as sns
from itertools import izip
import matplotlib.pyplot as plt
sns.set_style('white')
font = {'family' : 'normal',
        'sans-serif':'Arial',
        'size'   : 25}
rc('font', **font)

figure_path = '/stor/work/Lambowitz/cdw2854/bisufite_seq/figures'
figurename = figure_path + '/methyl_qc.pdf'


# In[2]:

def plot_insert_size(ax):
    datapath = '/stor/work/Lambowitz/cdw2854/bisufite_seq/insertSize'
    # BS-seq insert size table
    meth_df = pd.read_table(datapath+'/P13B_mix_S1.tsv')        .sort_values('isize')        .assign(data_type = 'Bisulfite-treated cfDNA')        .pipe(lambda d: d[['isize','percentage','data_type']])

    # normal cfDNA
    plasma_df = pd.read_csv('/stor/work/Lambowitz/cdw2854/plasmaDNA/insertSize/P1022_1_S3_umi2id.csv')        .groupby(['isize'], as_index=False)        .agg({'count':'sum'})        .assign(percentage = lambda d: np.true_divide(d['count'],d['count'].sum()))        .sort_values('isize')        .assign(data_type = 'Plasma cfDNA')        .pipe(lambda d: d[['isize','percentage','data_type']])


    df = pd.concat([plasma_df,meth_df],axis=0)        .sort_values('isize')        .query('isize > 25')
    # plot isize
    color = ["windows blue", "amber", "greyish", "faded green", "dusty purple"]
    pal = sns.xkcd_palette(color)
    for datatype, col in izip(df['data_type'].unique(),pal):
        d = df.query('data_type=="%s"' %(datatype))
        ax.plot(d['isize'],d['percentage'],color = col,label=datatype)
    ax.legend(title = ' ',loc='upper right')
    ax.set_xlabel('Fragment size')
    ax.set_ylabel('% fragments')
    ax.set_yticklabels(ax.get_yticks() *100)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


# In[3]:

def plot_methyl_dens(ax):
    # read methylation sites
    methyl_plasma = '/stor/work/Lambowitz/cdw2854/bisufite_seq/methyl_calling/P13B_mix_S1_CpG.meth.bedGraph'
    plasma_methyl_df = pd.read_table(methyl_plasma,
                  names = ['chrom','start','end','methyl_dens'],
                  skiprows=1)
    plasma_methyl_df.head()

    # plot methylation density
    with sns.plotting_context('paper', font_scale=1.5):
        p = sns.distplot(plasma_methyl_df['methyl_dens']*100,
                         hist=True, kde=True, bins=50, ax =ax )
    plt.yticks(p.get_yticks(), p.get_yticks()*100)
    p.set_xlabel('Methylation density')
    p.set_ylabel('Percentage')
    p.spines['top'].set_visible(False)
    p.spines['right'].set_visible(False)


# In[8]:

plt.figure(figsize=(10,5))
ax = plt.subplot(121)
plot_insert_size(ax)
ax = plt.subplot(122)
plot_methyl_dens(ax)
plt.savefig(figurename)
print 'Saved: ', figurename


# In[ ]:



