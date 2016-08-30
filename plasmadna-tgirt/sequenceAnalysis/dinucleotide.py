#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import sys
import os
import re
import pandas as pd
import glob
from pybedtools import BedTool, set_tempdir
from pybedtools.cbedtools import Interval
import seaborn as sns
from multiprocessing import Pool
from scipy.signal import medfilt
from itertools import product
sns.set_style('white')

def allDinucleotides():
    bases = ['A','C','T','G','N']
    dinucleotides_header = [''.join(combination) for combination in product(bases,repeat=2)]
    return dinucleotides_header

def makeDinucleotideMatrix(dinucleotides, window):
    length = window * 2
    dinucleotide_dict = {di : {i: 0 for i in range(length)} for di in dinucleotides} 
    return dinucleotide_dict

def extractDinucleotide(dinucleotide_count_dict, sequence):
    first_nucleotide = sequence[:-1] 
    second_nucleotide = sequence[1:] 
    dinucleotides = map(lambda x,y: x+y, first_nucleotide, second_nucleotide)
    for i,di in enumerate(dinucleotides):
        dinucleotide_count_dict[di][i] += 1
    return dinucleotide_count_dict

def makeDF(nucleotideDict, lenType, window):
    df = pd.DataFrame.from_dict(nucleotideDict) \
        .assign(rowsum = lambda d: d.sum(axis=1))\
        .assign(position = lambda d: d.index - window) \
        .pipe(pd.melt, id_vars=['position','rowsum'], var_name='dinucleotide', value_name = 'count')\
        .assign(fraction = lambda d: np.true_divide(d['count'],d['rowsum'])) \
        .assign(lenType = lenType)
    return df

def fragmentCenter(fragment,window):
    center = (long(fragment.end) + long(fragment.start)) / 2
    center_start = center - window
    center_end = center + window
    insert_size = fragment.score
    return Interval(chrom = fragment.chrom,
            start = center_start, 
            end = center_end,
            score = insert_size,
            strand = fragment.strand)

def parseBed(bed_file, ref_fasta, window, regular_chrom):
    dinucleotides_header = allDinucleotides()
    short_dinucleotides_count = makeDinucleotideMatrix(dinucleotides_header, window)
    long_dinucleotides_count = makeDinucleotideMatrix(dinucleotides_header, window)
    large_fragments = [165,170]
    small_fragments = [35,80]
    long_count, short_count = 0, 0
    for fragment in BedTool(bed_file)\
            .filter(lambda x: x.chrom in regular_chrom)\
            .each(fragmentCenter, window)\
            .nucleotide_content(fi=ref_fasta, s =True, seq =True):
        insert_size = int(fragment.score)
        sequence = fragment.fields[-1]
        if small_fragments[0] <= insert_size <= small_fragments[1]:
            short_dinucleotides_count = extractDinucleotide(short_dinucleotides_count, sequence)
            short_count += 1
        elif large_fragments[0] <= insert_size <= large_fragments[1]:
            long_dinucleotides_count = extractDinucleotide(long_dinucleotides_count, sequence)
            long_count += 1
    print 'Parsed %i long fragments and %i short fragments' %(long_count, short_count)
    iterable = zip([short_dinucleotides_count, long_dinucleotides_count],
                ['short (%i-%ibp)' %(small_fragments[0],small_fragments[1]),
                'long (%i-%ibp)' %(large_fragments[0],large_fragments[1])])
    dfs = [makeDF(nucleotideDict, lenType, window) for nucleotideDict, lenType in iterable]
    df = pd.concat(dfs)
    return df

def renameDinucleotide(di):
    ''' A function to rename dinucleotides '''
    AT = list(map(lambda x: ''.join(x), product(['A','T'],repeat=2)))
    CG = list(map(lambda x: ''.join(x), product(['G','C'],repeat=2)))
    if di in AT:
        new = 'AT|TA|TT|AA'
    elif di in CG:
        new = 'CG|GC|CC|GG'
    else:
        new = 'Nope'
    return new

def plotting(df, figurename):
    #start plotting
    df = df[df['lenType'].str.contains('long')]
    with sns.plotting_context('paper', font_scale=2.5, rc={"lines.linewidth":2}):
        p = sns.FacetGrid(data = df, hue = 'dinucleotide_type',
                col = 'samplename', size = 4, aspect = 1.6, col_wrap=3) 
        p.map(plt.plot, 'position', 'deviation')
        p.set_titles('{col_name}', fontweight='bold')
        p.set_ylabels('log2(Observer/Expected)')
        p.set_xlabels('Postion relative to middle of transcripts')
        p.set_xticklabels(rotation=80)
        p.add_legend()
        p.savefig(figurename)
    return 0

def medianFilter(dt):
    dt['deviation'] = dt['fraction'] / medfilt(dt['fraction'],201)
    return dt

def runFile(args):
    bed_file, ref_fasta, window, outputpath = args
    samplename = os.path.basename(bed_file).split('.')[0]
    outputPrefix = outputpath + '/' + samplename
    figurename = outputPrefix + '.pdf'
    tablename = outputPrefix + '.tsv'
    regular_chrom = map(str, np.arange(1,23))
    regular_chrom.extend(['X','Y'])
    print 'Start analyzing %s ' %samplename
    df = parseBed(bed_file, ref_fasta, window, regular_chrom)\
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
    bedFilePath = projectpath + '/rmdup_bed'
    outputpath = projectpath + '/nucleotidesAnaylsis/dinucleotides'
    makedir(outputpath)
    set_tempdir(outputpath)
    bedFiles = glob.glob(bedFilePath + '/*bed')
    bedFiles = filter(lambda x:re.search('DB|SRR|NT|RNase|PD',x), bedFiles)
    outputprefix = outputpath + '/dinucleotides'
    tablename = outputprefix + '.tsv'
    figurename = outputprefix + '.pdf'
    window = 400
    if not os.path.isfile(tablename):
        pool = Pool(24)
        dfs = pool.map(runFile,[(bed_file, reference, window, outputpath) for bed_file in bedFiles])
        pool.close()
        pool.join()
        df = pd.concat(dfs)
        df.to_csv(tablename,sep='\t', index=False)
    df = pd.read_csv(tablename,sep='\t')
    figurename = outputprefix + '.pdf'
    plotting(df, figurename)
    print 'Saved: %s.' %tablename

if __name__ == '__main__':
    main()
