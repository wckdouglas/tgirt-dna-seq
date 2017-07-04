#!/usr/bin/env python

import numpy as np
import pandas as pd
import os
import glob
import re
from itertools import izip
from functools import partial
from multiprocessing import Pool



def get_gene_table():
    gene_expression_table = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/rna.csv'
    gene_expression_table = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/rna_type.csv'
    ge = pd.read_csv(gene_expression_table)
    ge.columns = ['id','name','cells','TPM','unit','cat',
                  'tissue_type','cell_type','description','rname']
    ge = ge \
        .query('TPM > 0')\
        .groupby(['id'], as_index=False)\
        .agg({'TPM':'count'}) \
        .query('TPM >= 3') \
        .drop('TPM', axis=1)\
        .merge(ge, how='inner') 

    return ge

def read_file(ge, filename):
    low_p, hi_p = 190, 199
    return pd.read_table(filename) \
            .query('periodicity <= %i & periodicity >= %i' %(hi_p, low_p))\
            .groupby(['name','type','id'], as_index=False) \
            .agg({'intensity':'mean'})\
            .merge(ge,'inner')\
            .drop(['name','type','id'], axis = 1) \
            .groupby(['cells','tissue_type'])\
            .corr(method='pearson')\
            .reset_index() \
            .query('level_2 == "TPM"')\
            .drop(['level_2','TPM'],axis=1)\
            .rename(columns = {'intensity':'cor'})  \
            .assign(cor = lambda d: -(d['cor']))  \
            .assign(samplename= os.path.basename(filename).replace('.tsv','')) 


def tissue_cor(selected_files, ge):
    reading_input = partial(read_file, ge)
    p = Pool(24)
    df = p.map(reading_input, selected_files)
    p.close()
    p.join()
    df = pd.concat(df, axis=0)\
            .pipe(pd.pivot_table, index=['cells','tissue_type'],
                        columns = 'samplename',
                        values = 'cor')\
            .reset_index()
    return df

def main():
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    ge = get_gene_table()
    for fft in ['R','scipy']:
        for gene_region in ['tss','gene_body']:
            datapath = projectpath + '/periodicity_' + gene_region +'/'+fft
            print datapath
            files = glob.glob(datapath + '/*tsv')
            selected_files = filter(lambda x: not re.search('cor', x), files)

            pdf = tissue_cor(selected_files, ge)
            tablename = datapath + '/cor_table.tsv'
            pdf.to_csv(tablename,sep='\t',index=False)
            print 'Written %s' %tablename

if __name__ == '__main__':
    main()
