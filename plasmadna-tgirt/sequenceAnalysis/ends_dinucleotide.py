#!/usr/bin/env python

from matplotlib import use as mpl_use
mpl_use('Agg')
import seaborn as sns
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
import glob
from multiprocessing import Pool
from pybedtools import BedTool, set_tempdir
from pybedtools.cbedtools import Interval
from functools import partial
from itertools import izip, product
sns.set_style('white')

def baseCountDict():
    baseDict = {''.join(base) : 0 for base in product(list('ACTGN'), repeat=2)}
    return baseDict

def createCountMat(nucleotides_window):
    countMat = {i: baseCountDict() for i in np.arange(nucleotides_window)}
    return countMat

def makeDF(seq_count_matrix, end, nucleotides_half_window):
    df = pd.DataFrame.from_dict(seq_count_matrix, orient='index')
    df['total'] = df.sum(axis=1)
    df['index'] = np.arange(nucleotides_half_window*2) - nucleotides_half_window
    for col in df.columns:
        df[col] = np.true_divide(df[col], df['total'])
    df['End'] = end
    return df

def plotNucleotideFreq(df, figurename):
    df = pd.melt(df, id_vars = ['index','End','sample'],
            var_name='Base',
            value_name = 'fraction')
    plt.figure(figsize=(15, 20))
    with sns.plotting_context('paper', font_scale=1.5):
        p = sns.FacetGrid(data=df, row='End', hue = 'Base', col = 'sample')
    p.map(plt.plot,'index','fraction', alpha=0.7)
    p.set_xlabels('Position relative to fragment ends')
    p.set_ylabels('Fraction of bases')
    p.set_titles('{row_name}: {col_name}')
    p.add_legend()
    [ax.axvline(x=0, ymin=0, ymax = 1, alpha = 0.9, color = 'silver') for ax in p.fig.axes]
    p.savefig(figurename, dpi=500)
    return 0

def makeBedLineWindow(bed_record, nucleotides_half_window):
    return Interval(chrom = bed_record.chrom,
            start = long(bed_record.start) - nucleotides_half_window,
            end = long(bed_record.end) + nucleotides_half_window,
            strand = bed_record.strand)

def caluculateEndsNucleotide(sequence, seq_mat, window_size, iter_window):
    seq5 = sequence[:(window_size+1)]
    seq3 = sequence[-(window_size+1):]
    seq5_next = seq5[1:]
    seq3_next = seq3[1:]
    iterable = izip(seq5, seq3, seq5_next, seq3_next, iter_window)
    for base_5, base_3, base_5_next, base_3_next, i in iterable :
        seq_mat["5'"][i][base_5 + base_5_next] += 1
        seq_mat["3'"][i][base_3 + base_3_next] += 1
    return seq_mat

def runFile(nucleotides_half_window, ref_fasta, outputpath, regularChromosome, bedFile):
    print 'Running %s' %bedFile
    samplename = os.path.basename(bedFile).split('.')[0]
    figurename = outputpath + '/' + samplename + '.pdf'
    tablename = figurename.replace('pdf','tsv')
    window_size = nucleotides_half_window * 2
    iter_window = np.arange(window_size)
    seq_mat = {end : createCountMat(window_size) for end in ["5'","3'"]}
    seq_count = 0
    for seq in BedTool(bedFile)\
            .filter(lambda x: x.chrom in regularChromosome)\
            .filter(lambda frag: 0 < long(frag.end) - long(frag.start) < 300)\
            .each(makeBedLineWindow, nucleotides_half_window) \
            .nucleotide_content(fi=ref_fasta, seq=True, s= True):
        sequence = seq.fields[-1]
        seq_mat = caluculateEndsNucleotide(sequence, seq_mat, window_size, iter_window)
        seq_count += 1
        if seq_count % 10000000 == 0:
            print 'Parsed %i fragments from %s' %(seq_count, samplename)
    dfs = [makeDF(seq_mat[end], end, nucleotides_half_window) for end in seq_mat.keys()]
    df = pd.concat(dfs)
    df['sample'] = samplename
    plotNucleotideFreq(df, figurename)
    df.to_csv(tablename, sep='\t',index=False)
    print 'Saved: %s and %s.' %(figurename, tablename)
    return df

def makedir(directory):
    if not os.path.isdir(directory):
        os.makedirs(directory)
    return 0

def main():
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    ref_path = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA'
    ref_fasta = ref_path + '/genome_rDNA.fa'
    bedFilePath = projectpath + '/bedFiles'
    outputpath = projectpath + '/nucleotidesAnaylsis/ends_diucleotides'
    outputprefix = outputpath + '/end_dinucleotide'
    tablename = outputprefix + '.tsv'
    figurename = outputprefix + '.pdf'
    bedFiles = glob.glob(bedFilePath + '/*.bed')
    nucleotides_half_window = 20
    regularChromosome = np.arange(1,22)
    regularChromosome = np.append(regularChromosome,['X','Y'])
    makedir(outputpath)
    set_tempdir(outputpath)
    func = partial(runFile, nucleotides_half_window, ref_fasta, outputpath, regularChromosome)
    dfs = Pool(24).map(func,bedFiles)
    df = pd.concat(dfs)
    df.to_csv(tablename,sep='\t', index=False)
    df = pd.read_csv(tablename,sep='\t')
    plotNucleotideFreq(df, figurename)
    print 'Saved: %s.' %figurename

if __name__ == '__main__':
    main()
