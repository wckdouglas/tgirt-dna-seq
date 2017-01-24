#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
from scipy import stats
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
    ssDNA_seq = 'SRR2130051.tsv'
    TGIRT_seq1 = 'P1022_3_S5.tsv'
    TGIRT_seq2 =  'P1203-SQ2_S3.tsv'
    tablename = datapath + '/cor_periodicity_table.csv'
    table_list = [[ssDNA_seq, TGIRT_seq1],
                [TGIRT_seq2, TGIRT_seq1]]
    labels = ['ssDNA-seq vs TGIRT-seq','TGIRT-seq vs TGIRT-seq']
    df = map(make_comparison_table, table_list, labels)
    df = pd.concat(df, axis=0)
    df.to_csv(tablename,index=False)
    print 'Written %s' %tablename

if __name__ == '__main__':
    main()
