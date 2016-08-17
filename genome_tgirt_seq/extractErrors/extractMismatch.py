#!/usr/bin/env python

import pysam
import numpy as np
import pandas as pd
import glob
import re
import os
from multiprocessing import Pool
from pybedtools import BedTool
import matplotlib
matplotlib.use('Agg')
import  matplotlib.pyplot as plt
import seaborn as sns
discard = 3

def cigarToSeq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: SMMMMMMMMMMMMMMMMIIIIIMMMMMMMMMMSSSS
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['M','I','D'],dtype='string')) #no softclipped
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def read_dict():
    mismatch_dict = {}
    bases = list('ACTGN')
    for ref in bases:
        mismatch_dict[ref] = {}
        for read in bases:
            mismatch_dict[ref][read] = 0
    return mismatch_dict

def extract_mismatch(aln, mismatch_dict, samplename):
    read = aln.query_alignment_sequence.upper()
    quality = aln.query_alignment_qualities
    ref = aln.get_reference_sequence().upper()
    cigar = aln.cigarstring
    cigar_seq = cigarToSeq(cigar)
    assert len(cigar_seq.replace('D','')) == len(read), 'Wrong cigar conversion!'
    fixed_read = []
    fixed_qual = []
    offset = 0
    for i in np.arange(len(cigar_seq)):
        cb = cigar_seq[i]
        if cb == 'M':
            fixed_read.append(read[i - offset])
            fixed_qual.append(quality[i - offset])
        elif cb == 'D':
            fixed_read.append('N')
            offset += 1
    assert len(fixed_read) == len(ref), 'Wrong fix read sequence\n' + ''.join(cigar_seq)+'\n' +''.join(fixed_read) +'\n'+ ref
    seq_len = len(read)
    for i, read_base, ref_base, read_q in zip(range(seq_len),read, ref, quality):
        if read_q > 30 and discard < i < seq_len - discard :
            mismatch_dict[ref_base][read_base] += 1
    #        if read_base != ref_base:
    #            print '%s: %s' %(samplename,aln.query_name)
    return mismatch_dict

def parse_bam(args):
    bam_file = args
    samplename = os.path.basename(bam_file).split('.')[0]
    print 'Analyzing %s' %samplename
    with pysam.Samfile(bam_file, 'rb') as bam:
        mismatch_dict = read_dict()
        for aln in bam:
            if not aln.is_unmapped and not aln.is_secondary and 'N' not in aln.cigarstring:
                mismatch_dict = extract_mismatch(aln, mismatch_dict, samplename)
    mismatch_df = pd.DataFrame(mismatch_dict)\
	    .assign(read = lambda d: d.index) \
	    .pipe(pd.melt, id_vars='read', var_name='ref', value_name='count')\
	    .assign(samplename = samplename)
    return mismatch_df

def norm_df(df):
    df['count'] = np.true_divide(df['count'],np.sum(df['count']))
    return df

def plot_mismatch(df, figure_name):
    columns = 2 if len(np.unique(df['samplename'])) > 1 else 1
    df['mismatch'] = map(lambda x,y: x + '-to-' + y, df['ref'],df['read'])
    df = df\
	.groupby(['samplename'])\
	.apply(norm_df)\
	.reset_index()\
	.pipe(lambda x: x[(x['ref']!='N') &(x['read'] != 'N') & (x['read']!=x['ref'])]) \
	.fillna(0)
    with sns.plotting_context('paper',font_scale=1.2):
	p = sns.FacetGrid(data = df, col = 'samplename', col_wrap = columns)
    p.map(sns.barplot, 'mismatch','count')
    p.set_xticklabels(labels = np.unique(df['mismatch']),rotation=90)
    p.set_titles('{col_name}')
    p.savefig(figure_name)

def main():
    project_path = '/scratch/02727/cdw2854/jurkatCells'
    bam_path = project_path + '/bamFiles'
    mismatch_path = project_path + '/mismatches'
    reference_path = '/corral-repl/utexas/2013lambowitz/Ref/hg19/Sequence/WholeGenomeFasta'
    bam_files = glob.glob(bam_path + '/*.bam')
    new_table_name = mismatch_path + '/mismatch.tsv'
    figure_name = new_table_name.replace('.tsv','.pdf')
    dfs = Pool(24).map(parse_bam, bam_files)
    df = pd.concat(dfs,axis=0)
    df.to_csv(new_table_name,sep='\t',index=False)
    print 'Written %s ' %new_table_name
    df = pd.read_table(new_table_name, sep='\t')
    plot_mismatch(df, figure_name)
    print 'Written %s ' %figure_name
    return 0

if __name__ == '__main__':
    main()
