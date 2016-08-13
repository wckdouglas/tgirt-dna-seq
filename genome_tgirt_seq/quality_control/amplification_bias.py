#!/usr/bin/env python

from Bio.SeqIO.QualityIO import FastqGeneralIterator
import pysam
import gzip
import numpy
import seaborn as sns
import glob
import matplotlib.pyplot as plt
import pandas as pd
import os
from multiprocessing import Pool
sns.set_style('white')

def parse_fastq(fastq_file):
    id_array = []
    member = []
    with gzip.open(fastq_file,'rb') as fq:
        for id_line, seq, qual in FastqGeneralIterator(fq):
            line = id_line.split(' ')
            id_array.append(line[0])
            member.append(line[1])
    df = pd.DataFrame({'id':id_array,'member_count':member})
    return df

def parse_bam(bam_file):
    id_array = []
    isize = []
    with pysam.Samfile(bam_file,'rb') as bam:
        for aln in bam:
            if aln.is_read1 and not aln.is_supplementary \
                    and not aln.is_secondary and aln.is_proper_pair:
                isize.append(abs(aln.template_length))
                id_array.append(aln.query_name)
    df = pd.DataFrame({'id':id_array,'isize':isize})\
        .pipe(lambda x: x[x['isize']>10])
    return df

def make_df(args):
    samplename, bam_path, fastq_path, suffix = args
    fastq_file = fastq_path + '/' + samplename + suffix
    bam_file = bam_path + '/' + samplename + '_0.bam'
    fastq_df = parse_fastq(fastq_file)
    bam_df = parse_bam(bam_file)
    df = fastq_df.merge(bam_df,how='inner')
    df['name'] = samplename
    return df

def plot_amplification(df, figurename):
    with sns.plotting_context('paper',font_scale=1.3):
	ax = sns.kdeplot(data = df['isize'],data2=df['member_count'],shade=True, cmap='Reds')
    ax.set(xlabel = 'Insert Size',ylabel = 'Barcode Family Member')
    plt.savefig(figurename)
    print 'Written %s' %figurename

def main():
    project_path = '/scratch/02727/cdw2854/jurkatCells'
    fastq_path = project_path + '/splitted'
    bam_path = project_path + '/bamFiles'
    table_path = project_path + '/figures'
    table_name = table_path + '/amplification_table.tsv'
    figure_name = table_name.replace('.tsv','.pdf')
    suffix = '_R1_001.fastq.gz'
    fastqs = glob.glob(fastq_path + '/*' + suffix)
    samplenames = map(lambda x: os.path.basename(x).replace(suffix,''), fastqs)
    args = map(lambda x: (x,bam_path, fastq_path,suffix),samplenames)
    df = Pool(24).map(make_df, args)
    df = pd.concat(df,axis=0)
    df.to_csv(table_name, sep='\t', index=False)
    print 'written %s' %table_name
    df = pd.read_table(table_name,sep='\t')
    plot_amplification(df, figure_name)


if __name__ == '__main__':
    main()
