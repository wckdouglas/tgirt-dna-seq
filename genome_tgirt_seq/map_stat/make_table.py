#!/usr/bin/env python

import pandas as pd
import numpy as np

def raw_df():
    seq_read = '/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12/CUT_75BP/seq_count.txt'
    raw_read = pd.read_table(seq_read, sep=':', names=['filename','read_count']) \
        .assign(samplename = lambda d: d.filename.str.replace('_R1_001.fastq.gz','')) \
        .pipe(lambda d: d[d.samplename.str.contains('UMI|NEB|kh|nextera')]) \
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
        .pipe(lambda d: d[d.samplename.str.contains('UMI|nextera|NEB|kh')]) \
        .pipe(lambda d: d[d.samplename.str.contains('umi|75')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('.stats',''))  \
        .assign(samplename = lambda d: d.samplename.str.replace('75bp_','')) \
        .pipe(lambda d: d[d.samplename.str.contains('umi|nextera')]) \
        .assign(samplename = lambda d: d.samplename.str.replace('_umi2id','')) 
    return map_read

def make_str(x):
    return ['%.1f' %(f * 100) for f in x]

def add_comma(x):
    x = int(x)
    b = ''
    for i, c in enumerate(str(x)[::-1]):
        b += c
        if (i+1) % 3 == 0:
            b += ','
    return b[::-1].strip(',')



def merge_columns(x,y):
    return x.apply(int).apply(add_comma) + '\n(' + make_str(y) + '%)' 

def get_sample_id(x):
    sample_id = 0
    if 'Nextera' in x['prep']:
        sample_id = x['samplename'].split('_')[0]
    elif 'UMI' in x['samplename']:
        sample_id = x['samplename'].split('_')[2]
    elif 'NEB' in x['samplename']:
        sample_id = x['samplename'].split('_')[1]
    else:
        sample_id = x['samplename'].split('_')[2]
    return sample_id

def get_fragmentation(x):
    prep = ''
    if 'nextera' in x:
        prep = 'Transposon'
    elif 'UMI' in x:
        prep = 'Covaris'
    elif 'kh' in x:
        prep = 'Covaris'
    elif 'NEB' in x:
        prep = 'NEB fragmentase'
    return prep
    
def get_R1(x):
    prep = ''
    if 'nextera' in x:
        prep = 'N/A'
    elif 'UMI' in x:
        prep = 'UMI only'
    elif 'kh' in x:
        prep = 'UMI w/ constant region'
    elif 'NEB' in x:
        prep = 'UMI w/ constant region'
    return prep

def main():
    map_read = map_df()
    raw_read = raw_df()
    df = map_read.merge(raw_read,how='inner') \
        .assign(prep = lambda d: map(lambda x: 'Nextera-XT' if 'nextera' in x else 'TGIRT-seq', d.samplename)) \
        .assign(umi_rate = lambda d: np.true_divide(d.umi, d.raw)) \
        .assign(supplementary_rate = lambda d: np.true_divide(d.supplementary,d.mapped)) \
        .assign(fusion = lambda d: merge_columns(d.supplementary, d.supplementary_rate)) \
        .assign(trimmed_rate = lambda d: np.true_divide(d['trimmed reads'], d.umi)) \
        .assign(mapping_rate = lambda d: np.true_divide(d.mapped, d['trimmed reads']*2)) \
        .assign(concordant_rate = lambda d: np.true_divide(d['proper pair'], d['mapped'])) \
        .assign(duplicate_rate = lambda d: np.true_divide(d['duplicates'], d['mapped'])) \
        .assign(umi = lambda d: merge_columns(d.umi, d.umi_rate)) \
        .assign(trimmed = lambda d: merge_columns(d['trimmed reads'], d.trimmed_rate)) \
        .assign(mapped = lambda d: merge_columns(d.mapped, d.mapping_rate)) \
        .assign(pair = lambda d: merge_columns(d['proper pair'], d.concordant_rate)) \
        .assign(duplicates = lambda d: merge_columns(d['duplicates'], d.duplicate_rate)) \
        .assign(sample_id = lambda d: [get_sample_id(r) for i, r in d.iterrows()]) \
        .assign(fragmentation = lambda d: map(get_fragmentation, d.samplename)) \
        .assign(R1_primer = lambda d: map(get_R1, d.samplename)) \
        .assign(raw = lambda d: map(add_comma, d.raw)) \
        .pipe(lambda d: d[['prep','sample_id','fragmentation','R1_primer','raw',
                           'umi','trimmed','mapped','pair','duplicates', 'fusion']])  \
        .sort_values(['prep','fragmentation','R1_primer','sample_id']) \
        .sort_values('prep',ascending = False)
    df.columns = ['Method','Sample ID','Fragmentation','R1 primer','Raw reads', 'UMI > Q20', 'Trimmed reads', 
                  'Mapped reads','Concordant pair','Duplication reads','Chimeric reads']
    tablename = '/stor/work/Lambowitz/cdw2854/ecoli_genome/figures/map_summary.csv'
    df.to_csv(tablename,index=False)
    print 'Written %s' %(tablename)

if __name__ == '__main__':
    main()




