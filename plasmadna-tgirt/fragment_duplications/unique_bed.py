#!/usr/bin/env python

import numpy as np
from scipy.spatial.distance import hamming
import fileinput
from itertools import combinations
import pandas as pd
from functools import partial
from collections import defaultdict


def hamming_barcode(a, b):
    assert len(a) == len(b), 'Wrong barcode extraction'
    return hamming(list(a), list(b)) * len(a)


def merge_record(df):
    barcodes_dict = defaultdict(set)
    member_dict = defaultdict(str)
    for i in range(len(df)):
        d = df.iloc[i,:]
        if d.distance < 2:
            if d.a not in member_dict.keys() and d.b not in member_dict.keys():
                barcodes_dict[d.a].add(d.b)
                member_dict[d.b] = d.a

            elif d.a in member_dict.keys() and d.b not in barcodes_dict.keys():
                barcodes_dict[member_dict[d.a]].add(d.b)

            elif d.b in member_dict.keys() and d.a not in barcodes_dict.keys():
                barcodes_dict[member_dict[d.b]].add(d.a)
    duplicate_barcodes = []
    [duplicate_barcodes.extend(list(b)) for b in barcodes_dict.itervalues()]
    return duplicate_barcodes


def make_hamming_df(barcodes):
    comparison = combinations(set(barcodes),r=2)
    df = pd.DataFrame(list(comparison), columns = ['a','b']) \
        .assign(distance = lambda d: map(hamming_barcode, d.a, d.b )) 
    [barcodes.remove(r) for r in merge_record(df)]
    return barcodes


def make_line(chrom, start, end, strand, barcode):
    length = long(end) - long(start)
    template = '{chrom}\t{start}\t{end}\t{barcode}\t{length}\t{strand}' \
                .format(chrom = chrom,
                        start = start,
                        end = end,
                        barcode = barcode,
                        length = length,
                        strand= strand)
    return template


for line in fileinput.input():
    line = line.rstrip()
    fields = line.split('\t')
    barcodes = fields[-1]
    template = partial(make_line, fields[0], fields[1], fields[2], fields[3])
    if ',' not in barcodes:
        print template(barcodes)
    else:
        barcodes = barcodes.split(',') 
        barcodes = list(set(barcodes))
        if len(barcodes) == 1:
            print template(barcodes[0])
        else:
            barcodes = make_hamming_df(barcodes)
            for barcode in barcodes:
                print template(barcode)



            
        
