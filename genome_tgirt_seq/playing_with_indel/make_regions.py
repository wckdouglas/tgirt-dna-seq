#!/usr/env python

import pandas as pd
import os
import sys
import numpy as np
from itertools import izip

indel_cut_off = int(sys.argv[1])
ref_path = os.environ['REF'] +'/Ecoli'
ref_table = ref_path + '/k12_mg1655_repeat_index.tsv'
df = pd.read_table(ref_table) \
    .assign(indel_index = lambda d: d.negative_index + d.positive_index) \
    .query('indel_index >= %i ' %indel_cut_off) 

out_file_name = ref_path + '/k12_mg1655_high_indel.bed'
count = 0
with open(out_file_name, 'w') as out:
    iterator =  df.iterrows()
    for i, base in iterator:
        if base['negative_index'] == base['indel_index']:
            start = base['start']
            mononucleotide = base['fwd_base']
            indel_index = base['indel_index']
            taken_base = 1
        elif taken_base != indel_index and base['fwd_base'] == mononucleotide:
            taken_base += 1
        elif taken_base == indel_index:
            assert base['positive_index'] == indel_index and base['fwd_base'] == mononucleotide,'Wrong parsing'
            end = base['start']
            line = 'NC_000913.3\t%i\t%i\tIndel%i\t%i\t+\t%s' %(start, end, count , indel_index, mononucleotide)
            out.write(line + '\n')
            count += 1
        else:
            print base
print 'Finsihed writing %s' %out_file_name
