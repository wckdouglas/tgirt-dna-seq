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
from multiprocessing import Pool
from collections import defaultdict
import os

complement = string.maketrans('ACTGNactgn','TGACNtgacn')

def plot_dist(dist):
    ax=plt.subplot(111)
    sns.barplot(dist.keys(),dist.values())
    ax.set_xticklabels(dist.keys(),rotation=90)


def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def extract_interval(ref_fasta, seq_length_dict, insert_dist, base_dist,
                    outputprefix, fold, chrom):
    sys.stderr.write('Simulating from chrom: %s\n' %(chrom))
    fasta = Fasta(ref_fasta)
    seq_count = 0
    total_length = seq_length_dict[chrom]
    starts = np.linspace(1, total_length, total_length / 1e6)
    ends = starts[1:]
    filename = outputprefix + '_' + chrom + '.bed'
    outfile = open(filename,'w')
    for s, e in izip(starts, ends):
        sys.stderr.write('Parsing %s: %i-%i\n' %(chrom,s,e))
        sequence = fasta.get_seq(chrom, int(s), int(e))
        for i in xrange(len(sequence) - 2 + 1):
            di_nucleotide_5 = str(sequence[i:i+ 2])
            out, out_rev = 0, 0
            if 'N' not in di_nucleotide_5:
                for cov in xrange(fold):
                    strand_watson = random.binomial(1, p = 0.5)
                    strand_creek = random.binomial(1, p = 0.5)
                    if strand_watson == 0:
                        out = random.binomial(1, p = base_dist["5'"][di_nucleotide_5])
                    elif strand_creek == 0:
                        out_rev = random.binomial(1, p = base_dist["5'"][reverse_complement(di_nucleotide_5)])

                    if out > 0:
                        insert_size = insert_dist.rvs()
                        start_site = s + i - 1
                        end_site = int(start_site + insert_size -1)
                        di_nucleotide_3 = str(fasta.get_seq(chrom, end_site, end_site + 1))
                        if 'N' not in di_nucleotide_3:
                            if random.binomial(1, p = base_dist["3'"][di_nucleotide_3]) == 1:
                                outfile.write( '%s\t%i\t%i\tSeq_%s_%i\t%i\t+\n' %(chrom, start_site, end_site,
                                                            chrom, seq_count, insert_size))
                                seq_count += 1

                    if out_rev > 0:
                        insert_size = insert_dist.rvs()
                        end_site = s + i 
                        start_site = int(end_site - insert_size -1)
                        di_nucleotide_3 = str(fasta.get_seq(chrom, start_site, start_site + 1))
                        if 'N' not in di_nucleotide_3:
                            if random.binomial(1, p = base_dist["3'"][di_nucleotide_3]) == 1:
                                outfile.write('%s\t%i\t%i\tSeq_%s_%i\t%i\t-\n' %(chrom, start_site, end_site,
                                                            chrom, seq_count, insert_size))
                                seq_count += 1
    outfile.close()
    sys.stderr.write('Written: %i for chrom: %s\n' %(seq_count, chrom))
    return filename


def profile_to_distribution(insert_profile_table, base_profile_table):
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))

    base_df = pd.read_csv(base_profile_table)

    base_dist = defaultdict(dict)
    iterable =  izip(base_df.End.values, 
                    base_df.dinucleotide.values, 
                    base_df.fraction.values)
    for end, di, frac in iterable:
        base_dist[end][di] = float(frac)
    return insert_dist, base_dist


def parse_opt():
    parser = argparse.ArgumentParser(description='Randomly select region as simulated tgirt-seq reads')
    parser.add_argument('-i','--insert', required=True, help = 'Insertion profile, csv file with 2 column named as: {isize and count}')
    parser.add_argument('-b','--base', required=True, help = "base profile, csv file with 2 column named as: {A,C,T,G, index (with respec to 5'end)}")
    parser.add_argument('-r','--refFasta', required=True, help = "Reference fasta")
    parser.add_argument('-f','--fold', type=int, default = 10, help = "how many fold of genome to simulate (default: 10)")
    parser.add_argument('-o','--outprefix', required=True, help = "outputprefix")
    parser.add_argument('-t','--threads', type=int, default = 1, help = 'Threads to use')
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
    seq_ids = [seq_id for seq_id in seq_ids if 'K' not in seq_id and
                                            'gi' not in seq_id and
                                            'GL' not in seq_id and
                                            'MT' not in seq_id]
    seq_length_dict = {seq_id:len(fasta[seq_id]) for seq_id in seq_ids}
    total_nucleotide = sum(seq_length_dict.values())
    seq_count = 0

    per_chrom_simulation = partial(extract_interval, ref_fasta, seq_length_dict,
                                insert_dist, base_dist, args.outprefix, fold)
    p = Pool(args.threads)
    outfiles = p.map(per_chrom_simulation, seq_ids)
    p.close()
    p.join()
    all_files = ' '.join(outfiles)
    command = 'cat %s > %s.bed' %(all_files, args.outprefix)
    os.system(command)
    command = 'rm %s' %(all_files)
    os.system(command)
    print 'Merged'


if __name__ == '__main__':
    main()
