#!/usr/bin/env python

from matplotlib import use
use('Agg', warn=False)
import pysam
import re
import os
import numpy as np
from multiprocessing import Pool
from collections import defaultdict
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
sns.set_style('white')

def extract_cigar(cigar_string):
    cigar_num = re.findall(r'[0-9]+',cigar_string)
    cigar_char = re.findall(r'[MIDSH]',cigar_string)
    return np.array(cigar_num, dtype='int'), np.array(cigar_char)
   
    
def find_position(aln, cigar_num, cigar_char):
    if (aln.is_read1 and aln.is_reverse) or (aln.is_read2 and not aln.is_reverse):
        cigar_num = cigar_num[::-1]
        cigar_char = cigar_char[::-1]
    cigar_num = cigar_num.cumsum()
    del_positions = cigar_num[cigar_char=='D']
    insertion_positions = cigar_num[cigar_char=='I']
    return del_positions, insertion_positions


def extract_indel_position(aln, indel_pos_dict):
    cigar_num, cigar_char = extract_cigar(aln.cigarstring)
    del_positions, insertion_positions = find_position(aln, cigar_num, cigar_char)
    read = 'Read 1' if aln.is_read1 else 'Read 2'
    for pos in del_positions:
        indel_pos_dict[read]['deletion'][pos]+=1
    for pos in insertion_positions:
        indel_pos_dict[read]['insertion'][pos] += 1
    return indel_pos_dict
    

def construct_dataframe(indel_pos_dict):
    dfs = []
    for read, indel_dict in indel_pos_dict.iteritems():
        for indel, position_dict in indel_dict.iteritems():
            df = pd.DataFrame({'position': position_dict.keys(),
                          'count': position_dict.values()}) \
                .assign(indel_type = indel) \
                .assign(read_type = read)
            if len(df) > 0:
                dfs.append(df)
    if len(dfs) > 1:
        return pd.concat(dfs, axis = 0)

def evaluate_aln(aln):
    pass_aln = not aln.is_supplementary and not aln.is_secondary
    pass_indel = 'I' in aln.cigarstring or 'D' in aln.cigarstring
    return pass_aln and pass_indel
    

def run_bam(bam_file):
    filename = os.path.basename(bam_file).split('.')[0]
    print 'Running %s' %(filename)
    indel_pos_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    with pysam.Samfile(bam_file,'rb') as bam:
        for aln in bam:
            if not aln.is_unmapped and evaluate_aln(aln):
                indel_pos_dict = extract_indel_position(aln, indel_pos_dict)
    df = construct_dataframe(indel_pos_dict)
    if df is not None:
        df = df.assign(samplename = filename)
        return df
    
def normalize_count(df):
    df['norm_count'] = df['count']/df['count'].sum()
    return df
    

def plot_indel(indel_df, figurename):
    indel_df = indel_df \
        .pipe(lambda d: d[d.samplename.str.contains('nextera|kh|kq')]) \
        .pipe(lambda d: d[d.samplename.str.contains('nextera|cluster')]) \
        .assign(prep = lambda d: map(lambda x: 'Nextera XT 'if 'nextera' in x else 'TGIRT-seq', d.samplename))\
        .groupby(['prep','position','read_type','indel_type'])\
        .agg({'count':np.sum})\
        .reset_index()\
        .groupby(['prep','read_type','indel_type'])\
        .apply(normalize_count) \
        .reset_index()
        
    with sns.plotting_context('paper',font_scale=1.3):
        p = sns.FacetGrid(data = indel_df, hue = 'prep',
                          margin_titles=True, legend_out=False,
                          col = 'indel_type', row= 'read_type')
    p.map(plt.plot, 'position', 'count')
    [plt.setp(ax.texts, text="") for ax in p.axes.flat]
    p.set_titles(row_template='{row_name}', 
            col_template="{col_name}",
            fontweight='bold', size=20)
    p.set_xlabels('Read Position')
    p.set_ylabels('Number of InDel')
    p.add_legend(title=' ')
    p.savefig(figurename)
    return 0


def main():
    project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome'
    bam_path = project_path +  '/picard_results'
    indel_pos_table_path = project_path + '/indel_pos_tables'
    if not os.path.isdir(indel_pos_table_path):
        os.makedirs(indel_pos_table_path)
    tablename = indel_pos_table_path + '/indel_position.tsv'
    figurename = tablename.replace('.tsv','.pdf')
    bam_files = glob.glob(bam_path + '/*.MarkDuplicate.bam')
    p = Pool(12)
    #df = p.map(run_bam, bam_files)
    p.close()
    p.join()
    #df = filter(lambda x: x is not None,df)
    #df = pd.concat(df, axis = 0)
    #df.to_csv(tablename, index = False, sep='\t')
    print 'Written %s' %tablename
    df = pd.read_table(tablename)
    plot_indel(df, figurename)
    print 'Written %s' %figurename
    

if __name__ == '__main__':
    main()
