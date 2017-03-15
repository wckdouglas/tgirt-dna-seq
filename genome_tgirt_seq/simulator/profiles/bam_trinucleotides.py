#!/usr/bin/env python

from matplotlib import use
use('Agg')
import matplotlib.pyplot as plt
import pysam
import numpy as np
from collections import defaultdict
import seaborn as sns
import pandas as pd
import string
from itertools import izip
import sys

complement = string.maketrans('ACTGN','TGACN')
def reverse_complement(seq):
    return seq.translate(complement)[::-1]

def norm_data(d):
    d['trinucleotide_fraction'] = np.true_divide(d.trinucleotide_count, d.trinucleotide_count.sum())
    return d

def make_dataframe(nucleotide_dict, end):
    return pd.DataFrame.from_dict(nucleotide_dict[end], orient = 'index') \
        .reset_index() \
        .rename(columns = {'index':'trinucleotides',
                           0:'trinucleotide_count'})\
        .assign(read_end = end)  \
        .fillna(0)

def plot_ends(df, figurename):
    with sns.plotting_context('paper',font_scale = 1.2), \
            sns.color_palette("husl", 8):
        p = sns.FacetGrid(data = df, col = 'read_end',
                      hue = 'base', aspect = 1.5)
    p.map(plt.plot, 'positions','base_fraction')
    p.add_legend()
    p.set_titles('{col_name}')
    p.set_axis_labels('Positions','Fraction')
    p.savefig(figurename)
    print 'Written %s ' %figurename
    return 0


def extract_nucleotides(bam):
    end_nucleotide_dict = defaultdict(lambda : defaultdict(int))
    for count, aln in enumerate(bam):
        condition_1 = (not aln.is_unmapped and not aln.is_supplementary)
        condition_2 = (not aln.is_duplicate and aln.mapping_quality > 1)
        if condition_1 and condition_2:
            read = "5'" if aln.is_read1 else "3'"
            sequence = str(aln.query_alignment_sequence)
            sequence = sequence if not aln.is_reverse else reverse_complement(sequence)
            #sequence = sequence.translate(complement)[::-1] if aln.is_read2 else sequence
            trinucleotide = str(sequence[:3])
            end_nucleotide_dict[read][trinucleotide] += 1
        if count % 10000000 == 0 and count != 0:
            print 'Parsed %i alignments' %(count)
    return end_nucleotide_dict


def main():
    if len(sys.argv) != 3:
        sys.exit('[usage] python %s <bamfile> <outprefix>' %(sys.argv[0]))

    bam_file = sys.argv[1]
    outprefix = sys.argv[2]
    figurename = outprefix + '.pdf'
    tablename = outprefix + '.csv'
    with pysam.Samfile(bam_file,'rb') as bam:
        end_nucleotide_dict = extract_nucleotides(bam)
    df = pd.concat([make_dataframe(end_nucleotide_dict, end) for end in ["5'","3'"]])\
            .drop_duplicates()
    #plot_ends(df, figurename)
    df.to_csv(tablename, index=False)
    print 'Written %s' %tablename
    return 0



if __name__ == '__main__':
    main()
