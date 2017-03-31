import random
from pyfaidx import Fasta
from collections import defaultdict
import numpy as np
from itertools import product, izip
import string
import pandas as pd
from functools import partial
from scipy.stats import rv_discrete, bernoulli
import sys
from libc.stdlib cimport rand, RAND_MAX

cdef double RAND_SCALE = 1.0 / RAND_MAX
cpdef double next_rand():
    '''
    https://groups.google.com/forum/#!topic/cython-users/jc-3UK2Ffoc
    '''
    return rand() * RAND_SCALE


class bernoulli_generator:
    '''
    Fast bernoulli random variable generator
    '''
    def __init__(self, p):
        self.p = p

    def get_prob(self):
        out = 0 if next_rand() <= self.p else 1
        return out


cpdef str generate_line(chrom, start, end, seq_count, insert_size, strand):
    return '{chrom}\t{start_site}\t{end_site}\tSeq_{chrom}_{seq_count}\t{isize}\t{strand}'\
                .format(chrom = chrom, start_site = long(start), end_site = long(end),
                        seq_count = seq_count.value, isize = insert_size,
                        strand = strand)

complement = string.maketrans('ACTGNactgn','TGACNTGACN')
cpdef str reverse_complement(str sequence):
    return sequence.translate(complement)[::-1]



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


def kmer_generator(kmer_length):
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
        base_dist[end][kmer] = bernoulli_generator(p)

    # store distribution in hash table
    kmers_3 = kmer_generator(kmer_3)
    kmers_5 = kmer_generator(kmer_5)
    for _kmer in kmers_5:
        _kmer = ''.join(_kmer)
        if side == '3' or side == 'no':
            base_dist["5'"][_kmer] = bernoulli_generator(1)
        elif side == '5' or side == 'both':
            base_dist["5'"][_kmer] = base_dist["5'"].get(_kmer, bernoulli_generator(0))

    for _kmer in kmers_3:
        _kmer = ''.join(_kmer)
        if side == '3' or side == 'both':
            base_dist["3'"][_kmer] = base_dist["3'"].get(_kmer, bernoulli_generator(0))
        elif side == '5' or side == 'no':
            base_dist["3'"][_kmer] = bernoulli_generator(1)

    return base_dist


def len_profile(insert_profile_table):

    # make insert distribution
    insert_df = pd.read_csv(insert_profile_table)\
        .assign(px = lambda d: np.true_divide(d['count'].values,d['count'].values.sum()))
    insert_dist = rv_discrete(name='custm', values=(insert_df.isize, insert_df.px))
    return insert_dist

def extract_interval(side, ref_fasta, insert_profile_table, base_profile_table, kmer_5, kmer_3,
                    outprefix, fold, chrom, seq_count, iterable):

    cdef:
        int position, chrom_len
        int out_count = 0
        int strand, insert_size
        int out
        str line

    iternum, (start_chrom, end_chrom) = iterable
    max_kmer = max(kmer_5, kmer_3)
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
                    out = base_dist["5'"][k_nucleotide_5].get_prob()
                    if out == 1:
                        start_site = start_chrom + position - 1  #is for adjusting the 0-base python?
                        end_site = int(start_site + insert_size)
                        if end_site < (chrom_len - max_kmer):
                            k_nucleotide_3 = str(fasta.get_seq(chrom, end_site - kmer_3 + 1, end_site))
                            #assert len(k_nucleotide_3) == kmer, "Wrong extraction of + strand 3' kmer: " + k_nucleotide_3
                            if 'N' not in k_nucleotide_3 and base_dist["3'"][k_nucleotide_3].get_prob() == 1:
                                line = generate_line(chrom, start_site, end_site, seq_count, insert_size, '+')
                                seq_count.value += 1
                                outfile.write(line + '\n')
                else:
                    reverse_k_nucleotide_5  = reverse_complement(k_nucleotide_5)
                    out = base_dist["5'"][reverse_k_nucleotide_5].get_prob()
                    if out == 1:
                        end_site = start_chrom + position + 2 #reversed This is the start when the read is reversed
                        start_site = int(end_site - insert_size) #this is the end site
                        if start_site > 0:
                            k_nucleotide_3 = str(fasta.get_seq(chrom, start_site, start_site + kmer_3 - 1))
                            k_nucleotide_3 = reverse_complement(k_nucleotide_3)
                            #assert len(k_nucleotide_3) == kmer, "Wrong extraction of - strand 3' kmer: " + k_nucleotide_3
                            if 'N' not in k_nucleotide_3 and base_dist["3'"][k_nucleotide_3].get_prob() == 1:
                                line = generate_line(chrom, start_site - 1, end_site, seq_count, insert_size, '-')
                                seq_count.value += 1
                                outfile.write(line + '\n')
        if position % 1000000 == 0 and position != 0:
            print 'Parsed %i positions' %(position)
        if seq_count.value % 100000 == 0 and seq_count.value > 0:
            print 'Output %i reads' %(seq_count.value)
    outfile.close()
    return outfile_name
