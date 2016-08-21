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
from functools import partial
sns.set_style('white')

def parse_fastq(fastq_file):
    id_array = []
    member = []
    with gzip.open(fastq_file,'rb') as fq:
        for id_line, seq, qual in FastqGeneralIterator(fq):
            line = id_line.split('_')
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

def make_df(bam_path, fastq_path, suffix, samplename):
    fastq_file = fastq_path + '/' + samplename + suffix
    bam_file = bam_path + '/' + samplename + '_0.bam'
    fastq_df = parse_fastq(fastq_file)
    bam_df = parse_bam(bam_file)
    df = fastq_df.merge(bam_df,how='inner')
    df['name'] = samplename
    return df

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/jurkatCells'
    fastq_path = '/stor/work/Lambowitz/Data/NGS/JA16381/combined/splitted'
    bam_path = project_path + '/bamFiles'
    table_path = project_path + '/figures'
    table_name = table_path + '/amplification_table.tsv'
    suffix = '_R1_001.fastq.gz'
    fastqs = glob.glob(fastq_path + '/DB*' + suffix)
    samplenames = map(lambda x: os.path.basename(x).replace(suffix,''), fastqs)
    func = partial(make_df, bam_path, fastq_path, suffix)
    p = Pool(24)
    df = map(func, samplenames)
    p.close()
    p.join()
    df = pd.concat(df,axis=0)
    df.to_csv(table_name, sep='\t', index=False)
    print 'written %s' %table_name
    df = pd.read_table(table_name,sep='\t')

if __name__ == '__main__':
    main()
