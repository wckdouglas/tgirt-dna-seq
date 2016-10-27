#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
sns.set_style('white')

datapath = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/periodicity_tables'

def read_wps_table(tablename):
    samplename = 'ssDNA-seq' if 'SRR' in tablename else 'TGIRT-seq'
    df = pd.read_table(datapath + '/'+ tablename) \
        .pipe(lambda d: d[['name','periodicity','intensity','samplename']])  \
        .assign(prep = samplename)
    return df
    
tablename = ['SRR2130052.tsv', 'PD_merged.tsv']
df = map(read_wps_table, tablename)
df = pd.concat(df, axis = 0) \
    .pipe(pd.pivot_table, 
            columns = ['prep'], 
            index='name',
            values='periodicity') \
    .reset_index()\
    .pipe(lambda d: d[d['TGIRT-seq']!=0])\
    .pipe(lambda d: d[d['ssDNA-seq']!=0]) 
    
p = sns.FacetGrid(data = df)
p.map(sns.kdeplot, 'TGIRT-seq', 'ssDNA-seq', 
                    shade = True, cmap='Reds')
low_periodicity, high_periodicity = 150, 250
p.set(xlim = (low_periodicity, high_periodicity), 
    ylim=(low_periodicity, high_periodicity))
p.set(xlabel = 'Periodicity (TGIRT-seq)')
p.set(ylabel = 'Periodicity (ssDNA-seq)')
p.set(title='Periodicity in 10000bp window')
figurename = datapath+'/periodicity_cor.pdf'
p.savefig(figurename)
print 'Written:', figurename
