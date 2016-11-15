#!/usr/bin/env python

from matplotlib import use
use('Agg')
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
from multiprocessing import Pool, Manager
from collections import defaultdict
import os

complement = string.maketrans('ACTGNactgn','TGACNtgacn')

def plot_dist(dist, outprefix):
    d = pd.DataFrame.from_dict(dist) \
        .reset_index() \
        .rename(columns = {'index':'tri_nucleotides'})  \
        .pipe(pd.melt, id_vars = ['tri_nucleotides'],
            var_name = 'ends', value_name='fraction') \
        .assign(third = lambda d: map(lambda x: x[2],d.tri_nucleotides))\
        .assign(first = lambda d:map(lambda x: x[0],d.tri_nucleotides))

    with sns.plotting_context('paper'):
        p = sns.FacetGrid(data = d, row = 'ends',
                    sharex=False, sharey=False, aspect= 2, row_order = ["5'","3'"])
    p.map(sns.barplot, 'tri_nucleotides', 'fraction')
    [ax.set_xticklabels(ax.get_xticklabels(),rotation=90, fontsize = 7) for ax in p.fig.axes]
    figurename = outprefix + '.pdf'
    p.savefig(figurename)
    sys.stderr.write('Written %s\n' %figurename)


def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def extract_interval(side, ref_fasta, insert_profile_table, base_profile_table,
                    outprefix, fold, chrom, seq_count, iterable):
    iternum, (s, e) = iterable
    insert_dist, base_dist = profile_to_distribution_bayes(insert_profile_table, base_profile_table, side)
    #plot_dist(base_dist, outprefix)
    random.seed(iternum)
    fasta = Fasta(ref_fasta)
    chrom_len = len(fasta[chrom])
    sys.stderr.write('Parsing %s: %i-%i\n' %(chrom,s,e))
    sequence = fasta.get_seq(chrom, int(s), int(e))
    outfile_name = outprefix + '.' + str(iternum) + '.bed'
    outfile = open(outfile_name, 'w')
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
                        if end_site < chrom_len:
                            #tri_nucleotide_3 = str(fasta.get_seq(chrom, end_site - 2, end_site))
                            tri_nucleotide_3 = str(fasta.get_seq(chrom, end_site - 2, end_site-2))
                            if 'N' not in tri_nucleotide_3:
                                if random.binomial(1, p = base_dist["3'"][tri_nucleotide_3]) == 1:
                                    line = '%s\t%i\t%i\tSeq_%s_%i\t%i\t+\n' %(chrom, start_site, end_site,
                                                                chrom, seq_count.value, insert_size)
                                    outfile.write(line)
                                    seq_count.value += 1

                else:
                    out = random.binomial(1, p = base_dist["5'"][reverse_tri_nucleotide_5])
                    if out == 1:
                        insert_size = insert_dist.rvs()
                        end_site = s + i + 2 #reversed This is the start when the read is reversed
                        start_site = int(end_site - insert_size) #this is the end site
                        if start_site > 0:
                            #reverse_tri_nucleotide_3 = reverse_complement(str(fasta.get_seq(chrom, start_site, start_site + 2)))
                            reverse_tri_nucleotide_3 = reverse_complement(str(fasta.get_seq(chrom, start_site, start_site )))
                            if 'N' not in reverse_tri_nucleotide_3:
                                if random.binomial(1, p = base_dist["3'"][reverse_tri_nucleotide_3]) == 1:
                                    line = '%s\t%i\t%i\tSeq_%s_%i\t%i\t-\n' %(chrom, start_site - 1, end_site,
                                                                chrom, seq_count.value, insert_size)
                                    outfile.write(line)
                                    seq_count.value += 1
    outfile.close()
    return outfile_name


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



def conditional_probit(d):
    d['probit'] = d.trinucleotide_count/ d.trinucleotide_count.sum()
    return d

def merge_probit(df1, df2):
    return pd.merge(df1, df2, how ='inner')

def make_nucleotides(x,y):
    if x == "3'":
        return y[0]
    if x == "5'":
        return y

def profile_to_distribution_bayes(insert_profile_table, base_profile_table, side):
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))

    base_df = pd.read_csv(base_profile_table) \
        .pipe(lambda d: d[~d.trinucleotides.str.contains('N')]) \
        .groupby('read_end')\
        .apply(conditional_probit) \
        .assign(nucleotides = lambda d: map(make_nucleotides, d.read_end, d.trinucleotides)) \
        .groupby(['read_end','nucleotides'])\
        .agg({'probit':np.sum})\
        .reset_index()

    base_dist = defaultdict(lambda: defaultdict(float))
    for t,p,r in izip(base_df.nucleotides, base_df.probit, base_df.read_end):
        base_dist[r][t] = p

    return insert_dist, base_dist



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
    threads = args.threads
    outbed = args.outbed
    outprefix = outbed .split('.')[0]
    side = args.side

    fasta = Fasta(ref_fasta)
    seq_id = fasta.keys()[0]
    total_length = len(fasta[seq_id])

    starts = np.linspace(1, total_length, threads)
    ends = starts[1:]
    p = Pool(threads)
    seq_count = Manager().Value('i',0)
    per_site_simulation = partial(extract_interval, side, ref_fasta, insert_profile_table,
                                  base_profile_table, outprefix, fold, str(seq_id), seq_count)
    iterable = enumerate(zip(starts,ends))
    outfiles = p.map(per_site_simulation, iterable)
    p.close()
    p.join()
    all_files = ' '.join(outfiles)
    command = 'cat %s > %s.bed' %(all_files, outbed)
    os.system(command)
    command = 'rm %s' %(all_files)
    os.system(command)
    sys.stderr.write('Written %s\n' %outbed)

if __name__ == '__main__':
    main()
