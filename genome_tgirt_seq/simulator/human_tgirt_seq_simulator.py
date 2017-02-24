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


def extract_interval(ref_fasta, seq_length_dict, insert_profile_table, base_profile_table,
                    outputprefix, fold, side, iterable):
    (seed, chrom) = iterable
    insert_dist, base_dist = profile_to_distribution(insert_profile_table, base_profile_table, side)
    random.seed(seed)
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
        for i in xrange(len(sequence) - 3):
            tri_nucleotide_5 = str(sequence[i:i+ 3])
            reverse_tri_nucleotide_5  = reverse_complement(tri_nucleotide_5)
            if 'N' not in tri_nucleotide_5:
                for cov in xrange(fold):
                    strand = random.binomial(1, p = 0.5) # positive = 0, negative = 1
                    if strand == 0:
                        out = random.binomial(1, p = base_dist["5'"][tri_nucleotide_5])
                        if out == 1:
                            insert_size = insert_dist.rvs()
                            start_site = s + i  - 1  #is for adjusting the 0-base python?
                            end_site = int(start_site + insert_size)
                            if end_site < total_length:
                                tri_nucleotide_3 = str(fasta.get_seq(chrom, end_site - 2, end_site))
                                if 'N' not in tri_nucleotide_3:
                                    if random.binomial(1, p = base_dist["3'"][tri_nucleotide_3]) == 1:
                                        line = '%s\t%i\t%i\tSeq_%s_%i\t%i\t+\n' %(chrom, start_site, end_site,
                                                                chrom, seq_count, insert_size)
                                        outfile.write(line)
                                        seq_count += 1

                    else:
                        out = random.binomial(1, p = base_dist["5'"][reverse_tri_nucleotide_5])
                        if out == 1:
                            insert_size = insert_dist.rvs()
                            end_site = s + i + 2 #reversed This is the start when the read is reversed
                            start_site = int(end_site - insert_size) #this is the end site
                            if start_site > 0:
                                reverse_tri_nucleotide_3 = reverse_complement(str(fasta.get_seq(chrom, start_site, start_site + 2)))
                                if 'N' not in reverse_tri_nucleotide_3:
                                    if random.binomial(1, p = base_dist["3'"][reverse_tri_nucleotide_3]) == 1:
                                        line = '%s\t%i\t%i\tSeq_%s_%i\t%i\t-\n' %(chrom, start_site - 1, end_site,
                                                                chrom, seq_count, insert_size)
                                        outfile.write(line)
                                        seq_count += 1
    outfile.close()
    sys.stderr.write('Written: %i for chrom: %s\n' %(seq_count, chrom))
    return filename


def get_prob(base_df, side, pos, nuc, end):
    if side == '3':
        pos = 20 - pos + 1
    elif side == '5':
        pos = pos
    elif side == 'both':
        pos = pos if end == "5'" else 20 - pos + 1
    elif side == 'no':
        pos = pos if end == "3'" else 20 - pos + 1

    end_is_right = (base_df.end == end)
    pos_is_right = (base_df.pos == pos)
    nucleotide_is_right = (base_df.base == nuc)
    d = base_df[end_is_right & nucleotide_is_right & pos_is_right]
    return float(d.base_fraction)


def profile_to_distribution(insert_profile_table, base_profile_table, side):
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))

    base_df = pd.read_csv(base_profile_table)

    base_dist = defaultdict(lambda: defaultdict(float))

    extract_prob = partial(get_prob, base_df, side)
    for tri_nucleotides in  product('ACTG',repeat=3):
        tri_nucleotides = ''.join(tri_nucleotides)
        for end in ["5'","3'"]:
            base_fraction = [extract_prob(pos+1, nuc, end) for pos, nuc in enumerate(tri_nucleotides)]
            base_dist[end][tri_nucleotides] = np.prod(base_fraction)

    return insert_dist, base_dist


def parse_opt():
    parser = argparse.ArgumentParser(description='Randomly select region as simulated tgirt-seq reads')
    parser.add_argument('-i','--insert', required=True, help = 'Insertion profile, csv file with 2 column named as: {isize and count}')
    parser.add_argument('-b','--base', required=True, help = "base profile, csv file with 2 column named as: {A,C,T,G, index (with respec to 5'end)}")
    parser.add_argument('-r','--refFasta', required=True, help = "Reference fasta")
    parser.add_argument('-f','--fold', type=int, default = 10, help = "how many fold of genome to simulate (default: 10)")
    parser.add_argument('-o','--outprefix', required=True, help = "outputprefix")
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
    side = args.side


    fasta = Fasta(ref_fasta)
    seq_ids = fasta.keys()
    seq_ids = [seq_id for seq_id in seq_ids if 'K' not in seq_id and
                                            'gi' not in seq_id and
                                            'GL' not in seq_id and
                                            'M' not in seq_id and
                                            'MT' not in seq_id]
    seq_length_dict = {seq_id:len(fasta[seq_id]) for seq_id in seq_ids}
    total_nucleotide = sum(seq_length_dict.values())
    seq_count = 0

    per_chrom_simulation = partial(extract_interval, ref_fasta, seq_length_dict,
                                insert_profile_table, base_profile_table,
                                args.outprefix, fold, side)
    iterable = enumerate(seq_ids)
    p = Pool(args.threads)
    outfiles = p.map(per_chrom_simulation, iterable)
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
