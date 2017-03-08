#!/usr/bin/env python

import fileinput
from collections import defaultdict
import pandas as pd
import numpy as np
from functools import partial
import sys

if len(sys.argv) != 2:
    sys.exit('usage: python %s <outputprefix>' %sys.argv[0])

prefix = sys.argv[1]

def normalize_position(d):
    d['base_fraction'] = np.true_divide(d['base_count'], d['base_count'].sum())
    return d

def main():
    '''
    Fragment ends!!! not read ends for 3'!!
    '''
    positions = 20
    discard = 0
    end_nuc_dict = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
    len_dict = defaultdict(int)
    for sequence in sys.stdin:
        sequence = sequence.rstrip().upper()
        isize = len(sequence)
        if isize > positions:
            len_dict[str(isize)] += 1
            end_5 = sequence[:positions]
            end_3 = sequence[-positions:]
            for i in xrange(positions):
                end_nuc_dict["5'"][i][end_5[i]] += 1
                end_nuc_dict["3'"][i][end_3[i]] += 1
        else:
            discard += 1
    print 'Discated %i fragments' %(discard) 
    dfs = []
    for end, nucleotide_dict in end_nuc_dict.iteritems():
        for position, position_dict in nucleotide_dict.iteritems():
            d = pd.DataFrame({'base':position_dict.keys(),
                          'base_count':position_dict.values()}) \
                .assign(end = end) \
                .assign(position = position)
            dfs.append(d)

    pd.concat(dfs,axis = 0)\
        .groupby(['end','position'])\
        .apply(normalize_position) \
        .to_csv('%s_base_profile.csv' %prefix,index=False)

    pd.DataFrame({'isize':len_dict.keys(), 
                'count': len_dict.values()})\
        .assign(isize = lambda d: np.array(d.isize,dtype=np.int64)) \
        .query('isize < 600') \
        .to_csv('%s_len_profile.csv' %(prefix), index=False)
    print 'Written table'

if __name__ == '__main__':
    main()
