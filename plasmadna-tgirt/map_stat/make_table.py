#!/usr/bin/env python

import pandas as pd
import numpy as np


def add_comma(x):
    x = int(x)
    b = ''
    for i, c in enumerate(str(x)[::-1]):
        b += c
        if (i+1) % 3 == 0:
            b += ','
    return b[::-1].strip(',')


def raw_df():
    seq_read = '/stor/work/Lambowitz/cdw2854/plasmaDNA/raw_data/seq_count.tsv'
    raw_read = pd.read_table(seq_read, names=['filename','read_count']) \
        .assign(samplename = lambda d: d.filename.str.replace('_R1_001.fastq.gz','')) \
        .drop(['filename'],axis=1) \
        .assign(umi = lambda d: map(lambda x: 'umi' if 'umi2id' in x else 'raw', d.samplename)) \
        .assign(samplename = lambda d: d.samplename.str.replace('_umi2id','')) \
        .pipe(pd.pivot_table,index= 'samplename', columns='umi', values='read_count')  \
        .reset_index()\
        .assign(umi = lambda d: d.umi.fillna(d.raw))  \
        .query('raw > 0') \
        .pipe(lambda d: d[d.samplename.str.contains('P1022|P1016|P13|P1113|51$|52$')]) \
        .pipe(lambda d: d[~d.samplename.str.contains('clustered$')])
    return raw_read

def map_df():
    bam_read = '/stor/work/Lambowitz/cdw2854/plasmaDNA/bamFiles/stats/combined_stats.tsv'
    map_read = pd.read_table(bam_read) \
        .assign(samplename = lambda d: d.samplename.str.replace('.stats',''))  \
        .assign(samplename = lambda d: d.samplename.str.replace('75bp_','')) \
        .pipe(lambda d: d[d.samplename.str.contains('P1022|P1016|P13|P1113|51$|52$')]) \
        .pipe(lambda d: d[d.samplename.str.contains('umi|51$|52$')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('_umi2id',''))  \
        .pipe(lambda d: d[~d.samplename.str.contains('clustered$')])
    return map_read

def make_str(x):
    return ['%.1f' %(f * 100) for f in x]

def merge_columns(x,y):
    return x.apply(int).apply(add_comma) + '\n(' + make_str(y) + '%)' 

def make_prep(x):
    return map(lambda i: 'ssDNA-seq' if 'SRR' in i else 'TGIRT-seq', x)

def get_R1(x):
    R1 = ''
    if x.startswith('P1022'):
        R1 = 'UMI w/ constant regions'
    elif x.startswith('SRR'):
        R1 = 'N/A'
    else:
        R1 = 'UMI only'
    return R1


def main():
    map_read = map_df()
    raw_read = raw_df()
    df = map_read.merge(raw_read,how='inner') \
        .assign(umi_rate = lambda d: np.true_divide(d.umi, d.raw)) \
        .assign(trimmed_rate = lambda d: np.true_divide(d['trimmed reads'], d.umi)) \
        .assign(mapping_rate = lambda d: np.true_divide(d.mapped, d['trimmed reads']*2)) \
        .assign(concordant_rate = lambda d: np.true_divide(d['proper pair'], d['mapped'])) \
        .assign(duplicate_rate = lambda d: np.true_divide(d['duplicates'], d['mapped'])) \
        .assign(linked_read_rate = lambda d: np.true_divide(d['supplementary'],d['mapped'])) \
        .assign(umi = lambda d: merge_columns(d.umi, d.umi_rate)) \
        .assign(trimmed = lambda d: merge_columns(d['trimmed reads'], d.trimmed_rate)) \
        .assign(mapped = lambda d: merge_columns(d.mapped, d.mapping_rate)) \
        .assign(pair = lambda d: merge_columns(d['proper pair'], d.concordant_rate)) \
        .assign(duplicates = lambda d: merge_columns(d['duplicates'], d.duplicate_rate)) \
        .assign(translocation_read = lambda d: merge_columns(d['supplementary'],d.linked_read_rate))\
        .assign(R1 = lambda d: map(get_R1, d.samplename)) \
        .assign(prep = lambda d: make_prep(d.samplename)) \
        .assign(raw=lambda d: map(add_comma,d.raw))\
        .pipe(lambda d: d[['prep','samplename','R1','raw','umi','trimmed','mapped','pair','translocation_read']])\
        .sort_values(['prep','samplename','R1'])
    df.columns = ['Method','Sample ID','R1 primer','Raw reads', 'UMI > Q20', 'Trimmed reads', 
                  'Mapped reads','Concordant pair','Chimeric reads']
    tablename = '/stor/work/Lambowitz/cdw2854/plasmaDNA/figures/plasma_map_summary.csv'
    df.to_csv(tablename,index=False)
    print 'Written %s' %(tablename)

if __name__ == '__main__':
    main()




