#!/usr/bin/env python

import pandas as pd
import glob 

def read_data(filename):
    if 'celline' in filename:
        cat = 'cells'
    elif 'tissue' in filename:
        cat = 'tissues'
    return pd.read_csv(filename) \
            .assign(cat = cat)

def expression(file_path):
    files = glob.glob(file_path+'*.csv.zip')
    df = map(read_data, files)
    df = pd.concat(df,axis=0)
    return df

def cell_labels(file_path):
    link = file_path + '/labels.txt'
    return pd.read_table(link) \
        .rename(columns = {'Name':'tissues'})

def assign_cat(x, y):
    out = ''
    if pd.isnull(x):
        if y == 'tissues':
            out = 'Primary Tissue'
    else:
        out = x
    return out

file_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/'
labels = cell_labels(file_path)
tpm_table = expression(file_path)

cells = tpm_table[['cat','Sample']]\
        .drop_duplicates() \
        .rename(columns = {'Sample':'tissues'})\
        .merge(labels, on='tissues', how='left')\
        .assign(Category = lambda d: map(assign_cat,d.Category, d.cat)) \
        .query('Category != ""')  \
        .pipe(lambda d: d[~pd.isnull(d.Category)]) 
        
tpm_table= tpm_table \
        .rename(columns = {'Sample':'tissues'})\
        .merge(cells, how = 'inner', on=['tissues','cat']) 


out_file = file_path + '/rna_type.csv'
tpm_table.to_csv(out_file, index=False)
print 'Made %s' %out_file

