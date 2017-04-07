#!/usr/bin/env python

import numpy as np
import sys
import os
import re
import pandas as pd
import glob
from pybedtools import BedTool, set_tempdir
from pybedtools.cbedtools import Interval
from pyfaidx import Fasta
from multiprocessing import Pool
from scipy.signal import medfilt
from itertools import product
from functools import partial
import math
import pyximport
pyximport.install()
from sequence_analytic_tools import extract_dinucleotide

def allDinucleotides():
    bases = ['A','C','T','G','N']
    dinucleotides_header = [''.join(combination) for combination in product(bases,repeat=2)]
    return dinucleotides_header

def makeDinucleotideMatrix(dinucleotides, window):
    length = window * 2
    dinucleotide_dict = {di : {i: 0 for i in range(length)} for di in dinucleotides}
    return dinucleotide_dict



def makeDF(nucleotideDict, lenType, window):
    return pd.DataFrame.from_dict(nucleotideDict) \
        .assign(rowsum = lambda d: d.sum(axis=1))\
        .assign(position = lambda d: d.index - window) \
        .pipe(pd.melt, id_vars=['position','rowsum'], var_name='dinucleotide', value_name = 'count')\
        .assign(fraction = lambda d: np.true_divide(d['count'],d['rowsum'])) \
        .assign(lenType = lenType)

def fragmentCenter(fragment,window):
    center = (long(fragment.end) + long(fragment.start)) / 2
    center_start = center - window
    center_end = center + window
    return Interval(chrom = fragment.chrom,
            start = center_start,
            end = center_end,
            score = fragment.score,
            strand = fragment.strand)

def parseBed(bed_file, ref_fasta, window, regular_chrom):
    dinucleotides_header = allDinucleotides()
    long_dinucleotides_count = makeDinucleotideMatrix(dinucleotides_header, window)
    long_count, short_count = 0, 0
    for fragment in BedTool(bed_file)\
            .filter(lambda x: x.chrom in regular_chrom)\
            .each(fragmentCenter, window)\
            .nucleotide_content(fi=ref_fasta, s =True, seq =True):
        insert_size = int(fragment.score)
        sequence = fragment.fields[-1]
        if 167 == insert_size:
            long_dinucleotides_count = extract_dinucleotide(long_dinucleotides_count, sequence)
            long_count += 1
    print 'Parsed %i long fragments and %i short fragments' %(long_count, short_count)
    df = makeDF(long_dinucleotides_count, '167 bp', window)
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
    regular_chrom.extend(['X','Y'])
    print 'Start analyzing %s ' %samplename
    df = parseBed(bed_file, ref_fasta, window_size, regular_chrom)\
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
    plotting(df, figurename)
    print 'Saved %s' %figurename
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
    set_tempdir(outputpath)
    bedFiles = glob.glob(bedFilePath + '/*.bed')
    outputprefix = outputpath + '/pybed_dinucleotides'
    tablename = outputprefix + '.tsv'
    figurename = outputprefix + '.pdf'
    window_size = 400
    analyze_bam_func = partial(analyze_bam, reference, window_size, outputpath)
    if not os.path.isfile(tablename):
        pool = Pool(24)
        #dfs = pool.map(analyze_bam_func, bedFiles)
        dfs = map(analyze_bam_func, bedFiles)
        pool.close()
        pool.join()
        df = pd.concat(dfs)
        df.to_csv(tablename,sep='\t', index=False)
    df = pd.read_csv(tablename,sep='\t')
    figurename = outputprefix + '.pdf'
    print 'Saved: %s.' %tablename

if __name__ == '__main__':
    main()
