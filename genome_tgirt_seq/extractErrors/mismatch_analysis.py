#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import seaborn as sns
import glob
import os
from pybedtools import BedTool, set_tempdir
import argparse
sns.set_style('white')

def get_args():
    parser = argparse.ArgumentParser(description='Process nucleotide table and output mismatch counts')
    parser.add_argument('-i','--infile',required=True,help='Input table')
    parser.add_argument('-o','--outfile',required=True, help='Output table')
    parser.add_argument('-s','--snp', required=True, help='SNP bed file')
    parser.add_argument('-c','--cutoff', default = 0, help='Coverage cut off (defualt: 0)', type=int)
    args = parser.parse_args()
    return args

def infer_ref(record):
    record = record[list('ACTG')]
    if np.sum(record) * 0.8 < np.max(record):
        return np.argmax(record)
    else:
        return record.argmax()

def read_file(filename, out_file, snp_bed, coverage_cut_off):
    file_name = filename.split('/')[-1]
    samplename = '_'.join(file_name.split('_')[:-1])
    chrom = file_name.split('.')[0].split('_')[-1]
    columns = ['chrom','start','end','ref','coverage','strand']
    columns.extend(list('ACTGN'))
    file_size = os.path.getsize(filename)
    if file_size > 40:
        df = pd.read_table(filename) \
            .rename(columns={'refID':'chrom','refPos':'start','refBase':'ref'}) \
            .assign(end = lambda x: x['start']+2)\
            .assign(start = lambda x: x['start']+1)\
            .sort_values(by='start',axis=0)\
            .assign(strand = '+')\
            .pipe(lambda x: x[columns])\
            .pipe(lambda x: x[x['coverage']> coverage_cut_off]) 
        print 'Read %s' %file_name
        if df.shape[0] > 0:
            gdf =  BedTool()\
                .from_dataframe(df)\
                .intersect(b=snp_bed, v=True, sorted=True)\
                .to_dataframe(names = columns)\
                .assign(new_ref = lambda x: map(infer_ref, [x.iloc[i,:] for i in range(len(x))])) \
                .pipe(lambda x: x[x['new_ref']!='X']) \
                .drop(['start','end','chrom','strand'],axis=1) \
                .pipe(pd.melt,id_vars=['new_ref','ref','coverage'],
                    var_name='read',
                    value_name='count')\
                .assign(mismatch = lambda d: map(lambda x,y: x+' to '+y, d['new_ref'],d['read']))\
                .groupby(['mismatch','new_ref','read']) \
                .agg({'count':sum})\
                .reset_index()\
                .assign(samplename = samplename)\
                .to_csv(out_file, sep='\t',index=False)
            print 'Written %s' %out_file
            return out_file

def main(args):
    outfile = args.outfile
    infile = args.infile
    snp_bed = args.snp
    coverage_cut_off = args.cutoff
    result_path = os.path.dirname(outfile)
    set_tempdir(result_path)
    read_file(infile, outfile, snp_bed, coverage_cut_off)
    return 0

if __name__ == '__main__':
    args = get_args()
    main(args)
