#!/usr/bin/env python

from Bio import SeqIO
import pandas as pd
from multiprocessing import Pool
import os

def generate_index_array(sequence):
    current_base = 'X'
    base_count = 0

    repeat_index = []
    for i, base in enumerate(sequence):
        if base == current_base:
            base_count += 1
        else:
            current_base = base
            base_count = 0
        repeat_index.append(base_count)
        if i % 1000000 == 0:
            print 'Parsed %i position' %i
    return repeat_index


def main():
    ref_path = os.environ['REF'] + '/Ecoli'
    sequence = SeqIO.parse(ref_path + '/k12_mg1655.fa','fasta')\
        .next()\
        .seq

    positive_index, negative_index = Pool(2).map(generate_index_array, [sequence, sequence.reverse_complement()])
    df = pd.DataFrame({'start':xrange(len(sequence)),
                    'positive_index':positive_index,
                    'negative_index':negative_index[::-1],
                    'fwd_base':list(sequence)})
    df.to_csv(ref_path + '/k12_mg1655_repeat_index.tsv',sep='\t',index=False)


if __name__ == '__main__':
    main()
