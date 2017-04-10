
#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re
import pandas as pd
import glob
import seaborn as sns
from multiprocessing import Pool
from scipy.signal import medfilt
from itertools import product
from functools import partial
import math
import pyximport
pyximport.install()
from sequence_analytic_tools import parse_bed
sns.set_style('white')


def makeDF(nucleotideDict, lenType, window):
    return pd.DataFrame.from_dict(nucleotideDict) \
        .assign(rowsum = lambda d: d.sum(axis=1))\
        .assign(position = lambda d: d.index - window) \
        .pipe(pd.melt, id_vars=['position','rowsum'], var_name='dinucleotide', value_name = 'count')\
        .assign(fraction = lambda d: np.true_divide(d['count'],d['rowsum'])) \
        .assign(lenType = lenType)

def from_bed_to_df(bed_file, ref_fasta, window, regular_chrom):
    dinucleotide_count_dict = parse_bed(bed_file, ref_fasta, window, regular_chrom)
    df = makeDF(dinucleotide_count_dict, '167 bp', window)
    return df

AT = list(map(lambda x: ''.join(x), product(['A','T'],repeat=2)))
CG = list(map(lambda x: ''.join(x), product(['G','C'],repeat=2)))
AT_label = '|'.join(AT)
CG_label = '|'.join(CG)
def renameDinucleotide(di):
    ''' A function to rename dinucleotides '''
    if di in AT:
        label = AT_label
    elif di in CG:
        label = CG_label
    else:
        label = 'Nope'
    return label

def medianFilter(dt):
    dt['adjusted_signal'] = dt['fraction'] - medfilt(dt['fraction'],101)
    return dt

def analyze_bam(ref_fasta, window_size, outputpath, bed_file):
    samplename = os.path.basename(bed_file).split('.')[0]
    outputPrefix = outputpath + '/' + samplename
    figurename = outputPrefix + '.pdf'
    tablename = outputPrefix + '.tsv'
    regular_chrom = map(str, np.arange(1,23))
    print 'Start analyzing %s ' %samplename
    df = from_bed_to_df(bed_file, ref_fasta, window_size, regular_chrom)\
        .assign(dinucleotide_type = lambda d: map(renameDinucleotide, d['dinucleotide'])) \
        .pipe(lambda d: d[d['dinucleotide_type'] != 'Nope'])\
        .groupby(['dinucleotide_type','position','lenType'])\
        .agg({'fraction':np.sum})\
        .reset_index()\
        .groupby(['lenType','dinucleotide_type'])\
        .apply(medianFilter) \
        .reset_index() \
        .assign(samplename = samplename)
    df.to_csv(tablename, sep='\t',index=False)
    print 'Saved %s' %tablename
    return df

def makedir(directory):
    if not os.path.isdir(directory):
        os.system('mkdir -p %s' %directory)
    return 0

def main():
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    referencepath = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA'
    reference = referencepath + '/genome_rDNA.fa'
    bedFilePath = projectpath + '/bedFiles'
    outputpath = projectpath + '/nucleotidesAnaylsis/dinucleotides'
    makedir(outputpath)
    bedFiles = glob.glob(bedFilePath + '/*.bed')
    bedFiles = filter(lambda x: re.search('51_rmdup|52_rmdup|umi2id_unique|mix_unique',x), bedFiles)
    bedFiles = filter(lambda x: not re.search('SQ',x), bedFiles)
    bedFiles = filter(lambda x: re.search('mix',x), bedFiles)
    outputprefix = outputpath + '/dinucleotides'
    tablename = outputprefix + '.tsv'
    figurename = outputprefix + '.pdf'
    window_size = 400
    analyze_bam_func = partial(analyze_bam, reference, window_size, outputpath)
    pool = Pool(12)
    dfs = pool.map(analyze_bam_func, bedFiles)
    #dfs = map(analyze_bam_func, bedFiles)
    pool.close()
    pool.join()
    df = pd.concat(dfs)
    df.to_csv(tablename,sep='\t', index=False)
    print 'Saved: %s.' %tablename

if __name__ == '__main__':
    main()
