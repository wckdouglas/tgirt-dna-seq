#!/usr/bin/env python

from scipy.stats import rv_discrete, bernoulli
from scipy import random
import numpy as np
from pyfaidx import Fasta
from itertools import product, izip, cycle
from functools import partial
import string
import pandas as pd
import sys
import argparse
from multiprocessing import  Manager, Pool, Lock
from collections import defaultdict
import os


def generate_line(chrom, start, end, seq_count, insert_size, strand):
    return '{chrom}\t{start_site}\t{end_site}\tSeq_{chrom}_{seq_count}\t{isize}\t{strand}'\
                .format(chrom = chrom, start_site = long(start), end_site = long(end),
                        seq_count = seq_count.value, isize = insert_size,
                        strand = strand)

complement = string.maketrans('ACTGNactgn','TGACNTGACN')
def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def extract_interval(side, ref_fasta, insert_profile_table, base_profile_table, kmer_5, kmer_3,
                    outprefix, fold, chrom, seq_count, iterable):
    iternum, (start_chrom, end_chrom) = iterable
    max_kmer = max(kmer_5, kmer_3)
    out_count = 0
    insert_dist = len_profile(insert_profile_table)
    base_dist = base_profile(base_profile_table, side, kmer_5, kmer_3)
    print 'Constructed base distribution'
    #plot_dist(base_dist, outprefix)
    random.seed(iternum)
    fasta = Fasta(ref_fasta)
    chrom_len = len(fasta[chrom])
    sys.stderr.write('Parsing %s: %i-%i\n' %(chrom,start_chrom,end_chrom))
    sequence = fasta.get_seq(chrom, int(start_chrom), int(end_chrom))
    outfile_name = outprefix + '.' + str(iternum) + '.bed'
    outfile = open(outfile_name, 'w')
    for position in xrange(len(sequence) - max_kmer):
        k_nucleotide_5 = str(sequence[position: position+ kmer_5])
        #assert len(k_nucleotide_5) == kmer, "Wrong extraction of 5' kmer: " + k_nucleotide_5
        if 'N' not in k_nucleotide_5:
            strands = bernoulli.rvs(p = 0.5, size = fold)
            insert_sizes = insert_dist.rvs(size = fold)
            for strand, insert_size in zip(strands, insert_sizes):
                if strand == 1:
                    out = base_dist["5'"][k_nucleotide_5].next()
                    if out == 1:
                        start_site = start_chrom + position - 1  #is for adjusting the 0-base python?
                        end_site = int(start_site + insert_size)
                        if end_site < (chrom_len - max_kmer):
                            k_nucleotide_3 = str(fasta.get_seq(chrom, end_site - kmer_3 + 1, end_site))
                            #assert len(k_nucleotide_3) == kmer, "Wrong extraction of + strand 3' kmer: " + k_nucleotide_3
                            if 'N' not in k_nucleotide_3 and base_dist["3'"][k_nucleotide_3].next() == 1:
                                line = generate_line(chrom, start_site, end_site, seq_count, insert_size, '+')
                                seq_count.value += 1
                                outfile.write(line + '\n')
                else:
                    reverse_k_nucleotide_5  = reverse_complement(k_nucleotide_5)
                    out = base_dist["5'"][reverse_k_nucleotide_5].next()
                    if out == 1:
                        end_site = start_chrom + position + 2 #reversed This is the start when the read is reversed
                        start_site = int(end_site - insert_size) #this is the end site
                        if start_site > 0:
                            k_nucleotide_3 = str(fasta.get_seq(chrom, start_site, start_site + kmer_3 - 1))
                            k_nucleotide_3 = reverse_complement(k_nucleotide_3)
                            #assert len(k_nucleotide_3) == kmer, "Wrong extraction of - strand 3' kmer: " + k_nucleotide_3
                            if 'N' not in k_nucleotide_3 and base_dist["3'"][k_nucleotide_3].next() == 1:
                                line = generate_line(chrom, start_site - 1, end_site, seq_count, insert_size, '-')
                                seq_count.value += 1
                                outfile.write(line + '\n')
        if position % 1000000 == 0 and position != 0:
            print 'Parsed %i positions' %(position)
        if seq_count.value % 100000 == 0 and seq_count.value > 0:
            print 'Output %i reads' %(seq_count.value)
    outfile.close()
    return outfile_name

def define_scaling_factor(float_p):
    '''
    upscale all probability but limiting the highest prob < 1
    to save iteration for simulation
    '''
    return 10**(abs(int(np.log10(float_p))))

def extract_kmer(kmer_5, kmer_3, end, kmer):
    '''
    Trim kmer
    '''
    return kmer[:kmer_5] if end == "5'" else kmer[:kmer_3]


def prob_generator(kmer_length):
    '''
    making possible kmer
    '''
    kmers = product('ACTG',repeat=kmer_length)
    return kmers


def base_profile(base_profile_table, side, kmer_5, kmer_3):
    n = 1000
    # make kmer according to 5' or 3'
    get_kmer = partial(extract_kmer, kmer_5, kmer_3)
    base_df = pd.read_csv(base_profile_table) \
        .assign(kmer = lambda d: map(get_kmer, d['end'], d['kmer'])) \
        .groupby(['end','kmer'], as_index=False ) \
        .agg({'kmer_fraction':np.sum})

    fwd_max_p = base_df[base_df.end=="5'"]['kmer_fraction'].max()
    rvs_max_p = base_df[base_df.end=="3'"]['kmer_fraction'].max()
    fwd_scaling_factor = define_scaling_factor(fwd_max_p)
    rvs_scaling_factor = define_scaling_factor(rvs_max_p)

    # scale up kmers and make distribution
    base_dist = defaultdict(lambda: defaultdict(float))
    for i, row in base_df.iterrows():
        end = row['end']
        kmer = row['kmer']
        scaling_factor = fwd_scaling_factor if end == "5'" else rvs_scaling_factor
        p = row['kmer_fraction'] * scaling_factor
        base_dist[end][kmer] = cycle(bernoulli(p = p).rvs(n))

    # store distribution in hash table
    kmers_3 = prob_generator(kmer_3)
    kmers_5 = prob_generator(kmer_5)
    zero_generator = cycle([0])
    uniform_generator = cycle([1])
    for _kmer in kmers_5:
        _kmer = ''.join(_kmer)
        if side == '3' or side == 'no':
            base_dist["5'"][_kmer] = uniform_generator
        elif side == '5' or side == 'both':
            base_dist["5'"][_kmer] = base_dist["5'"].get(_kmer, zero_generator)

    for _kmer in kmers_3:
        _kmer = ''.join(_kmer)
        if side == '3' or side == 'both':
            base_dist["3'"][_kmer] = base_dist["3'"].get(_kmer, zero_generator)
        elif side == '5' or side == 'no':
            base_dist["3'"][_kmer] = uniform_generator

    return base_dist


def len_profile(insert_profile_table):

    # make insert distribution
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))
    return insert_dist


def parse_opt():
    parser = argparse.ArgumentParser(description='Randomly select region as simulated tgirt-seq reads')
    parser.add_argument('-i','--insert', required=True, help = 'Insertion profile, csv file with 2 column named as: {isize and count}')
    parser.add_argument('-b','--base', required=True, help = "base profile, csv file with 2 column named as: {A,C,T,G, index (with respec to 5'end)}")
    parser.add_argument('-r','--refFasta', required=True, help = "Reference fasta")
    parser.add_argument('-f','--fold', type=int, default = 10, help = "how many fold of genome to simulate (default: 10)")
    parser.add_argument('-o','--outbed', required=True, help = "outputprefix")
    parser.add_argument('-t','--threads', type=int, default = 1, help = 'Threads to use')
    parser.add_argument('-s','--side', choices = ['both','3','5','no'], required = True, help = 'which read to simulate with bias')
    args = parser.parse_args()
    return args


def main():
    args = parse_opt()
    insert_profile_table = args.insert
    base_profile_table = args.base
    ref_fasta = args.refFasta
    fold = args.fold
    threads = int(args.threads) + 1
    outbed = args.outbed
    outprefix = outbed .split('.')[0]
    side = args.side
    kmer_5, kmer_3 = 3, 12

    fasta = Fasta(ref_fasta)
    seq_id = fasta.keys()[0]
    total_length = len(fasta[seq_id])

    starts = np.linspace(1, total_length, threads)
    ends = starts[1:]
    p = Pool(threads)
    seq_count = Manager().Value('i',0)
#    base_dist = base_profile(base_profile_table, side, kmer_5, kmer_3)
#    print 'Constructed base distribution'
    per_site_simulation = partial(extract_interval, side, ref_fasta, insert_profile_table,
                                  base_profile_table, kmer_5, kmer_3, outprefix, fold, str(seq_id), seq_count)
    iterable = enumerate(zip(starts,ends))
    outfiles = p.map(per_site_simulation, iterable)
    #outfiles = map(per_site_simulation, iterable)
    p.close()
    p.join()
    with open(outbed + '.bed','w') as out:
        for f in outfiles:
            [out.write(line) for line in open(f,'r')]
    map(os.remove, outfiles)
    sys.stderr.write('Written %s\n' %outbed)

if __name__ == '__main__':
    main()
