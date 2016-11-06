#!/usr/bin

import numpy as np
from pybedtools import BedTool, set_tempdir
import pandas as pd
from functools import partial
from multiprocessing import Pool
from collections import defaultdict
from itertools import izip, imap
import os

def extract_dinucleotide(dinucleotide_counter, first_dinucleotide, second_dinucleotide,
                        end_dinucleotide, sequence):
    first_nucleotide = sequence[:-1]
    second_nucleotide = sequence[1:]
    seq_length = len(first_nucleotide)
    dinucleotides = imap(lambda x,y: x+y, first_nucleotide, second_nucleotide)
    for i, di in enumerate(dinucleotides):
        if i == 0:
            first_dinucleotide[di] += 1
        if i == 1:
            second_dinucleotide[di] += 1
        elif i == seq_length:
            end_dinucleotide[di] += 1

        dinucleotide_counter[di] += 1

    return dinucleotide_counter, first_dinucleotide,\
            second_dinucleotide, end_dinucleotide


def analyze_file(datapath, regular_chrom, ref_fasta, filename):
    print 'Running %s' %filename
    bed_name = datapath + '/' + filename
    dinucleotide_counter = defaultdict(int)
    first_dinucleotide = defaultdict(int)
    second_dinucleotide = defaultdict(int)
    end_dinucleotide = defaultdict(int)
    iterator =  BedTool(bed_name)\
        .filter(lambda f: f.chrom in regular_chrom) \
        .nucleotide_content(fi=ref_fasta, s =True, seq =True)
    print 'Make bed'

    for line_count, fragment in enumerate(iterator):
        sequence = fragment.fields[-1]
        if 'N' not in sequence:
            dinucleotide_counter, \
            first_dinucleotide, \
            second_dinucleotide, \
            end_dinucleotide = extract_dinucleotide(dinucleotide_counter,
                                                    first_dinucleotide,
                                                    second_dinucleotide,
                                                    end_dinucleotide,
                                                    sequence)
        if line_count % 1000000 == 0:
            print 'Parsed %i liens ' %(line_count)

    end_df = pd.DataFrame.from_dict(end_dinucleotide, orient='index') \
        .reset_index() \
        .rename(columns = {'index':'dinucleotide',0:'count'}) \
        .assign(what = 'end')
    second_df = pd.DataFrame.from_dict(second_dinucleotide, orient='index') \
        .reset_index() \
        .rename(columns = {'index':'dinucleotide',0:'count'}) \
        .assign(what = 'second')
    first_df = pd.DataFrame.from_dict(first_dinucleotide, orient='index') \
        .reset_index() \
        .rename(columns = {'index':'dinucleotide',0:'count'})  \
        .assign(what = 'first')
    dinucleotides_df = pd.DataFrame.from_dict(dinucleotide_counter, orient='index') \
        .reset_index() \
        .rename(columns = {'index':'dinucleotide',0:'count'}) \
        .assign(what = 'content')
    tablename = 'TGIRT.csv' if filename.startswith('P1203') else 'ssDNA-seq.csv'
    tablename = 'profiles/' + tablename
    df = pd.concat([first_df,second_df,end_df, dinucleotides_df])
    df.to_csv(tablename, index=False)
    print 'Written ', tablename

def main():
    bed_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA/bedFiles'
    set_tempdir(bed_path)
    ref_fasta = os.environ['REF'] + '/GRCh38/hg38_rDNA/genome_rDNA.fa'
    filenames = ['P1203-SQ2_S3.bed','SRR2130052.bed']
    regular_chrom = map(str, np.arange(1,23))
    regular_chrom.extend(['X','Y'])
    func = partial(analyze_file, bed_path, regular_chrom, ref_fasta)
    p = Pool(12)
    p.map(func, filenames)
    p.close()
    p.join()


if __name__ == '__main__':
    main()
