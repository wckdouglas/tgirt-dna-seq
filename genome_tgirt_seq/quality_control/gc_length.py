#!/usr/bin/env python

from __future__ import print_function
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from functools import partial
from multiprocessing import Pool
from pybedtools import BedTool, set_tempdir
from collections import defaultdict, Counter
import sys
from glob import glob
import os
sns.set_style('white')

def extractGClength(fragment):
    fields = fragment.fields
    return '_'.join([fields[-1],fields[7]])

def parseBed(resultpath, ref_fasta, bed_file):
    samplename = bed_file.split('/')[-1].split('.')[0]
    print('Running %s' %samplename, file=sys.stderr)
    df_name = resultpath + '/' + samplename + '.tsv'
    gc_len_dict = defaultdict(int)
    iterator = BedTool(bed_file) \
            .nucleotide_content(fi=ref_fasta,s=True, stream=True) 
    mat = defaultdict(int)
    for fragment in iterator:
        mat[extractGClength(fragment)]+=1
    df = pd.DataFrame(mat.items(), columns = ['seqlength_gc','count'])
    df = df['seqlength_gc'].str.split('_',expand = True) \
        .rename(columns={0:'seq_length',1:'gc'}) \
        .assign(count = df['count'])\
        .assign(samplename = samplename)
    print('saved %s' %df_name, file=sys.stderr)
    df.to_csv(df_name,sep='\t', index=False)
    return df

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bed_path = project_path + '/rmdup_bed'
    result_path = project_path + '/gc_length'
    ref_fasta = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA/genome_rDNA.fa'
    set_tempdir(result_path)
    if not os.path.isdir(result_path):
        os.makedirs(result_path)
    bed_files = glob(bed_path + '/*bed')
    df_name = result_path + '/gc_length.tsv'
    figurename = df_name.replace('.tsv','png')
    func = partial(parseBed, result_path, ref_fasta)
    p = Pool(12)
    dfs = p.map(func, bed_files)
    p.close()
    p.join()
    df = pd.concat(dfs, axis=0) 
    df.to_csv(df_name,sep='\t',index=False)
    print('Saved %s' %df_name, file=sys.stderr)

if __name__ == '__main__':
    main()
