#!/usr/bin/env python

import pandas as pd
import numpy as np

def raw_df():
    seq_read = '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/CUT_75BP/seq_count.txt'
    raw_read = pd.read_table(seq_read, sep=':', names=['filename','read_count']) \
        .assign(samplename = lambda d: d.filename.str.replace('_R1_001.fastq.gz','')) \
        .pipe(lambda d: d[d.samplename.str.contains('UMI|nextera')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('75bp_','')) \
        .drop(['filename'],axis=1) \
        .assign(umi = lambda d: map(lambda x: 'umi' if 'umi2id' in x else 'raw', d.samplename)) \
        .assign(samplename = lambda d: d.samplename.str.replace('_umi2id','')) \
        .pipe(pd.pivot_table,index= 'samplename', columns='umi', values='read_count')  \
        .reset_index()\
        .assign(umi = lambda d: d.umi.fillna(d.raw))  
    return raw_read

def map_df():
    bam_read = '/stor/work/Lambowitz/cdw2854/ecoli_genome/bamFiles/mark_duplicate/stats/combined_stats.tsv'
    map_read = pd.read_table(bam_read) \
        .pipe(lambda d: d[d.samplename.str.contains('UMI|nextera')]) \
        .pipe(lambda d: d[d.samplename.str.contains('umi|75')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('.stats',''))  \
        .assign(samplename = lambda d: d.samplename.str.replace('75bp_','')) \
        .pipe(lambda d: d[d.samplename.str.contains('umi|nextera')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('_umi2id','')) 
    return map_read

def make_str(x):
    return ['%.1f' %(f * 100) for f in x]

def merge_columns(x,y):
    return x.apply(str) + ' (' + make_str(y) + '%)' 

def get_sample_id(x):
    sample_id = ''
    if 'Nextera' in x['prep']:
        sample_id = x['samplename'].split('_')[0]
    elif 'UMI' in x['samplename']:
        sample_id = x['samplename'].split('_')[2]
    return sample_id

def main():
    map_read = map_df()
    raw_read = raw_df()
    df = map_read.merge(raw_read,how='inner') \
        .assign(prep = lambda d: map(lambda x: 'Nextera-XT' if 'nextera' in x else 'TGIRT-seq', d.samplename)) \
        .sort_values('prep') \
        .assign(umi_rate = lambda d: np.true_divide(d.umi, d.raw)) \
        .assign(trimmed_rate = lambda d: np.true_divide(d['trimmed reads']/2, d.umi)) \
        .assign(mapping_rate = lambda d: np.true_divide(d.mapped, d['trimmed reads'])) \
            .assign(concordant_rate = lambda d: np.true_divide(d['proper pair'], d['mapped'])) \
        .assign(duplicate_rate = lambda d: np.true_divide(d['duplicates'], d['mapped'])) \
        .assign(umi = lambda d: merge_columns(d.umi, d.umi_rate)) \
        .assign(trimmed = lambda d: merge_columns(d['trimmed reads'], d.trimmed_rate)) \
        .assign(mapped = lambda d: merge_columns(d.mapped, d.mapping_rate)) \
        .assign(pair = lambda d: merge_columns(d['proper pair'], d.concordant_rate)) \
        .assign(duplicates = lambda d: merge_columns(d['duplicates'], d.duplicate_rate)) \
        .assign(sample_id = lambda d: [get_sample_id(r) for i, r in d.iterrows()]) \
        .pipe(lambda d: d[['sample_id','prep','raw','umi','trimmed','mapped','pair','duplicates']]) 
    tablename = '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures/map_summary.csv'
    df.to_csv(tablename, index=False)
    print 'Written %s' %(tablename)

if __name__ == '__main__':
    main()




