#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
sns.set_style('white')

datapath = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/periodicity_tables'

def rename_sample(name):
    if 'SRR' in name:
        return 'ssDNA-seq'
    elif 'SQ1' in name:
        return 'TGIRT-seq 1'
    elif 'SQ2' in name:
        return 'TGIRT-seq 2'
    else:
        return 'TGIRT-seq merge'

def read_wps_table(tablename):
    samplename = rename_sample(tablename)
    print 'read %s' %samplename
    df = pd.read_table(datapath + '/'+ tablename) \
        .pipe(lambda d: d[['name','periodicity','intensity','samplename']])  \
        .rename(columns = {'periodicity':samplename})\
        .drop(['samplename','intensity'],axis = 1)
    return df
    

def merging(df1, df2):
    return pd.merge(df1, df2, how = 'inner', on = 'name')
    
def make_comparison_table(table_list, label):
    df = map(read_wps_table, table_list)
    df = reduce(merging, df) \
        .assign(minimum = lambda d: np.min(d.iloc[:,1:],axis = 1)) \
        .query('minimum > 0') \
        .drop('minimum',axis=1) 
    df.columns = ['name','sample1','sample2']
    df['label'] = label    
    return df

def main():
    table_list = [['SRR2130051.tsv', 'P1203_SQ2_S3_clustered.tsv'],
                ['P1203-SQ1_S2.tsv', 'P1203-SQ2_S3.tsv']]
    labels = ['ssDNA-seq vs TGIRT-seq','TGIRT-seq vs TGIRT-seq'] 
    df = map(make_comparison_table, table_list, labels)
    df = pd.concat(df, axis=0)
        
    low_periodicity, high_periodicity = 150, 250
    with sns.plotting_context('paper', font_scale=1.2):
        p = sns.FacetGrid(data = df, col = 'label')
    p.map(sns.kdeplot, 'sample1', 'sample2', 
                    shade = True, cmap='viridis')
    p.set(xlim=(low_periodicity, high_periodicity), 
        ylim=(low_periodicity, high_periodicity))
    p.set_axis_labels('Periodicity','Periodicity')
    p.set_titles('{col_name}')
    p.add_legend() 
    p.fig.suptitle('Nucleosome Spacing\n(10kb windows)',
        fontsize=14, fontweight='bold')
    p.fig.subplots_adjust(top=0.75)
    figurename = datapath + '/cor_periodicity.pdf'
    p.savefig(figurename)
    print 'Written %s' %figurename
    
if __name__ == '__main__':
    main()
