#!/usr/bin/env python

import numpy as np
from functools import partial
import string
import sys
import argparse
from multiprocessing import  Manager, Pool, Lock
from pyfaidx import Fasta
import os
import random
import pyximport
pyximport.install()
from simulator import extract_interval



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
