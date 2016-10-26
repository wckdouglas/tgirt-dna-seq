#!/usr/bin/env python

from scipy.stats import rv_discrete
from scipy import random
import numpy as np
from pyfaidx import Fasta
import seaborn as sns
import matplotlib.pyplot as plt
from itertools import product, izip
from functools import partial
import string
import pandas as pd
import sys
import argparse

complement = string.maketrans('ACTGNactgn','TGACNtgacn')

def plot_dist(dist):
    ax=plt.subplot(111)
    sns.barplot(dist.keys(),dist.values())
    ax.set_xticklabels(dist.keys(),rotation=90)


def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def extract_interval(fasta, seq_length_dict, insert_dist, base_dist, fold, seq_count, chrom):
    total_length = seq_length_dict[chrom]
    starts = np.linspace(1, total_length, total_length / 1e5)
    ends = starts[1:]
    for s, e in izip(starts, ends):
        sys.stderr.write('Parsing %s: %i-%i\n' %(chrom,s,e))
        sequence = fasta.get_seq(chrom, int(s), int(e))
        for i in xrange(len(sequence)-2):
            tri_nucleotide = str(sequence[i:i+3])
            if 'N' not in tri_nucleotide:
                out = sum(random.binomial(fold, base_dist[tri_nucleotide]))
                out_rev = sum(random.binomial(fold, base_dist[reverse_complement(tri_nucleotide)]))

                if out > 0:
                    insert_size = insert_dist.rvs()
                    start_site = s + i
                    end_site = start_site + insert_size
                    for y in xrange(out):
                        print '%s\t%i\t%i\tSeq_%i\t0\t+' %(chrom, start_site, end_site, seq_count)
                        seq_count += 1

                if out_rev > 0:
                    insert_size = insert_dist.rvs()
                    end_site = s + i + 2
                    start_site = end_site - insert_size
                    for y in xrange(out):
                        print '%s\t%i\t%i\tSeq_%i\t0\t-' %(chrom, start_site, end_site, seq_count)
                        seq_count += 1
    return seq_count


def profile_to_distribution(insert_profile_table, base_profile_table):
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))

    base_df = pd.read_csv(base_profile_table)

    possible_trinucleotide = product(list('ACTG'),repeat=3)
    base_dist = {''.join(tri): np.prod([base_df[b].values[i] for i, b in enumerate(tri)]) for tri in possible_trinucleotide}
    #plot_dist(base_dist)

    return insert_dist, base_dist


def parse_opt():
    parser = argparse.ArgumentParser(description='Randomly select region as simulated tgirt-seq reads')
    parser.add_argument('-i','--insert', required=True, help = 'Insertion profile, csv file with 2 column named as: {isize and count}')
    parser.add_argument('-b','--base', required=True, help = "base profile, csv file with 2 column named as: {A,C,T,G, index (with respec to 5'end)}")
    parser.add_argument('-r','--refFasta', required=True, help = "Reference fasta")
    parser.add_argument('-f','--fold', required=True, help = "how many fold of genome to simulate")
    args = parser.parse_args()
    return args


def main():
    args = parse_opt()
    insert_profile_table = args.insert
    base_profile_table = args.base
    ref_fasta = args.refFasta
    fold = args.fold
    insert_dist, base_dist = profile_to_distribution(insert_profile_table, base_profile_table)

    fasta = Fasta(ref_fasta)
    seq_ids = fasta.keys()
    seq_length_dict = {seq_id:len(fasta[seq_id]) for seq_id in seq_ids}
    total_nucleotide = sum(seq_length_dict.values())
    seq_count = 0

    per_chrom_simulation = partial(extract_interval, fasta, seq_length_dict, insert_dist, base_dist, np.arange(fold))
    for chrom in seq_ids:
        seq_count = per_chrom_simulation(seq_count, chrom)


if __name__ == '__main__':
    main()
