#!/usr/bin/env python

from scipy.stats import rv_discrete, bernoulli
from scipy import random
import numpy as np
from pyfaidx import Fasta
from itertools import product, izip
from functools import partial
import string
import pandas as pd
import sys
import argparse
from multiprocessing import Pool, Manager
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


def extract_interval(side, ref_fasta, insert_profile_table, base_profile_table,
                    outprefix, fold, chrom, seq_count, iterable):
    iternum, (start_chrom, end_chrom) = iterable
    kmer = 3
    out_count = 0
    insert_dist, base_dist = profile_to_distribution(insert_profile_table, base_profile_table, side, kmer)
    #plot_dist(base_dist, outprefix)
    random.seed(iternum)
    fasta = Fasta(ref_fasta)
    chrom_len = len(fasta[chrom])
    sys.stderr.write('Parsing %s: %i-%i\n' %(chrom,start_chrom,end_chrom))
    sequence = fasta.get_seq(chrom, int(start_chrom), int(end_chrom))
    outfile_name = outprefix + '.' + str(iternum) + '.bed'
    outfile = open(outfile_name, 'w')
    for position in xrange(len(sequence) - kmer):
        tri_nucleotide_5 = str(sequence[position: position+ kmer])
        #assert len(tri_nucleotide_5) == kmer, "Wrong extraction of 5' kmer: " + tri_nucleotide_5
        if 'N' not in tri_nucleotide_5:
            reverse_tri_nucleotide_5  = reverse_complement(tri_nucleotide_5)
            strands = bernoulli.rvs(p = 0.5, size = fold)
            insert_sizes = insert_dist.rvs(size = fold)
            for strand, insert_size in zip(strands, insert_sizes):
                if strand == 1:
                    out = base_dist["5'"][tri_nucleotide_5].rvs()
                    if out == 1:
                        start_site = start_chrom + position - 1  #is for adjusting the 0-base python?
                        end_site = int(start_site + insert_size)
                        if end_site < chrom_len - kmer:
                            tri_nucleotide_3 = str(fasta.get_seq(chrom, end_site - kmer + 1, end_site))
                            #assert len(tri_nucleotide_3) == kmer, "Wrong extraction of + strand 3' kmer: " + tri_nucleotide_3
                            if 'N' not in tri_nucleotide_3 and base_dist["3'"][tri_nucleotide_3].rvs() == 1:
                                line = generate_line(chrom, start_site, end_site, seq_count, insert_size, '+')
                                seq_count.value += 1
                                outfile.write(line + '\n')
                else:
                    out = base_dist["5'"][reverse_tri_nucleotide_5].rvs()
                    if out == 1:
                        end_site = start_chrom + position + 2 #reversed This is the start when the read is reversed
                        start_site = int(end_site - insert_size) #this is the end site
                        if start_site > 0:
                            tri_nucleotide_3 = reverse_complement(str(fasta.get_seq(chrom, start_site, start_site + kmer -1)))
                            #assert len(tri_nucleotide_3) == kmer, "Wrong extraction of - strand 3' kmer: " + tri_nucleotide_3
                            if 'N' not in tri_nucleotide_3 and base_dist["3'"][tri_nucleotide_3].rvs() == 1:
                                line = generate_line(chrom, start_site, end_site, seq_count, insert_size, '-')
                                seq_count.value += 1
                                outfile.write(line + '\n')
        if position % 1000000 == 0 and position != 0:
            print 'Parsed %i positions' %(position)
        if seq_count.value % 100000 == 0:
            print 'Output %i reads' %(seq_count.value)
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

def profile_to_distribution(insert_profile_table, base_profile_table, side, k):
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))

    base_df = pd.read_csv(base_profile_table)

    base_dist = defaultdict(lambda: defaultdict(float))

    extract_prob = partial(get_prob, base_df, side)
    for kmer in  product('ACTG',repeat=k):
        kmer = ''.join(kmer)
        for end in ["5'","3'"]:
            base_fraction = [extract_prob(pos+1, nuc, end) for pos, nuc in enumerate(kmer)]
            base_dist[end][kmer] = bernoulli(p = np.prod(base_fraction))

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
    #outfiles = map(per_site_simulation, iterable)
    p.close()
    p.join()
    all_files = ' '.join(outfiles)
    command = 'cat %s > %s.bed' %(all_files, outbed)
    os.system(command)
    #command = 'rm %s' %(all_files)
    #os.system(command)
    sys.stderr.write('Written %s\n' %outbed)

if __name__ == '__main__':
    main()
