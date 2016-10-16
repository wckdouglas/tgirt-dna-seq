#!/usr/bin/env python

import pandas as pd
import numpy as np
import glob

def retrieveCount(tablename):
    sample = tablename.split('/')[-1].split('.')[0]
    df = pd.read_table(tablename, 
                        names=['gene_id','read_length','count','unused']) \
        .drop(['unused','read_length'], axis = 1) \
        .pipe(lambda d: d[d['count']>0])  \
        .assign(samplename = sample)
    print 'Read %s' %sample
    return df

def normalizeCount(df):
    df['normalizedCount'] = np.true_divide(df['count'],np.sum(df['count']))
    return df

def getGeneTable():
    gene_table = '/stor/work/Lambowitz/ref/GRCh38/Bed_for_counts_only/genes.Info'
    tx_table = '/stor/work/Lambowitz/ref/human_transcriptome/tx.info'
    gene_df = pd.read_table(gene_table) 
    gene_df = pd.read_table(tx_table) \
        .rename(columns = {'gene':'GeneID'}) \
        .merge(gene_df, on = ['GeneID','type'], how='left') 
    return gene_df

def getHistone(name):
    return 'Histone Genes' if 'HIST' in name else 'MT-RNR' if ('MT-RNR' in name or 'MTRNR' in name) else 'Other genes'

def main():
    table_path = '/stor/work/Lambowitz/cdw2854/target-seq/removed_clipped'
    tables = glob.glob(table_path + '/*errorFree.tsv')
    new_table = table_path + '/count.tsv'
    df = map(retrieveCount, tables)
    df = pd.concat(df,axis = 0)\
        .groupby(['samplename','gene_id']) \
        .agg({'count':np.sum})\
        .reset_index() \
        .groupby(['samplename']) \
        .apply(normalizeCount) \
        .reset_index()

    df = getGeneTable() \
        .rename(columns = {'id':'gene_id'}) \
        .merge(df, on = ['gene_id'], how = 'right')  \
        .drop(['gene_id','index'],axis=1) \
        .fillna('name') \
        .assign(histone = lambda d: map(getHistone, d['name']))
    print 'Written %s' %new_table
    df.to_csv(new_table, sep='\t',index=False)

if __name__ == '__main__':
    main()
