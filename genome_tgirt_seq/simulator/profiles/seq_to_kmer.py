#!/usr/bin/env python

import fileinput
from collections import defaultdict
import pandas as pd
import numpy as np
from functools import partial
import sys

if len(sys.argv) != 3:
    sys.exit('usage: python %s <outputprefix> <k>' %sys.argv[0])

prefix = sys.argv[1]
k = int(sys.argv[2])
if k < 1 or k > 20:
    sys.exit('usage: python %s <outputprefix> <k>' %sys.argv[0])

def count_to_fraction(d):
    return d

def make_df(nuc_dict, key):
    return pd.DataFrame({'kmer':nuc_dict.keys(),
                         'kmer_count':nuc_dict.values()}) \
        .assign(end = key) \
        .assign(kmer_fraction = lambda d: np.true_divide(d.kmer_count,d.kmer_count.sum()))

def main():
    '''
    Fragment ends!!! not read ends for 3'!!
    '''
    end_nuc_dict = defaultdict(lambda: defaultdict(int))
    len_dict = defaultdict(int)
    for sequence in sys.stdin:
        sequence = sequence.rstrip().upper()
        isize = len(sequence)
        len_dict[str(isize)] += 1
        end_5_kmer = sequence[:k]
        end_3_kmer = sequence[-k:]
        end_nuc_dict["5'"][end_5_kmer] += 1
        end_nuc_dict["3'"][end_3_kmer] += 1
    
    dfs = [make_df(end_nuc_dict[key], key) for key in end_nuc_dict.iterkeys()]
    pd.concat(dfs,axis = 0)\
        .to_csv('%s_base_profile.csv' %prefix,index=False)

    pd.DataFrame({'isize':len_dict.keys(), 
                'count': len_dict.values()})\
        .assign(isize = lambda d: np.array(d.isize,dtype=np.int64)) \
        .query('isize < 600') \
        .to_csv('%s_len_profile.csv' %(prefix), index=False)
    print 'Written table'

if __name__ == '__main__':
    main()
