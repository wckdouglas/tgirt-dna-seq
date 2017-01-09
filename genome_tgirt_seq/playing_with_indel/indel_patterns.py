#!/usr/bin/env python

from matplotlib import use
use('Agg', warn=False)
import pysam
import re
import os
import sys
import numpy as np
from multiprocessing import Pool
from collections import defaultdict
import glob
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from itertools import imap, izip
from dmisc import sam_action
import string
sns.set_style('white')

def extract_cigar(cigar_string):
    cigar_num = re.findall(r'[0-9]+',cigar_string)
    cigar_char = re.findall(r'[MID]',cigar_string) # remove soft-clipped and hard-clip to match with reference seq
    cigar_string = ''
    for n, c in izip(cigar_num, cigar_char):
        cigar_string += int(n) * c
    return cigar_string


def normalize_sequence(cigar, md, seq):
    indel_indicator, corrected_seq = '',''
    for c, m, s in izip(cigar,md,seq):
        if m == 'D' and c == 'D':
            indel_indicator += 'D'
            corrected_seq += s
        elif c == 'I':
            indel_indicator += 'I'
            corrected_seq += 'I'
            corrected_seq += s
        else:
            corrected_seq += s
            indel_indicator += '-'
    return indel_indicator, corrected_seq


complement = string.maketrans('ACTGNactgn','TGACNTGACN')
def reverse_complement(sequence):
    return sequence.translate(complement)[::-1]


def adding_indel(indel_pattern_dict, indel_label, corrected_seq, positions):
#    indel_code = 'N'
#    if indel_label == 'Insertion':
#        indel_code = 'I' 
#    if indel_label == 'Deletion':
#        indel_code = 'D'
    seq_len = len(corrected_seq)
    limit = seq_len - 3
    rolling = 1
    for i in positions:
        if i < limit:
            indel_pattern = corrected_seq[i+1:i+4]
#            if indel_code in indel_pattern:
            if bool(re.search('[ID]',indel_pattern)):
                rolling += 1
            else: 
                indel_pattern_dict[indel_label][indel_pattern] += rolling
                rolling = 1
    return indel_pattern_dict

   
def analyze_alignment(aln, indel_pattern_dict):
    cigar_string = extract_cigar(aln.cigarstring)
    MD_seq = sam_action.MDToSeq(aln.get_tag('MD'))
    indel_indicator, corrected_seq = normalize_sequence(cigar_string, 
                                                        MD_seq, 
                                                        aln.get_reference_sequence())
    #print indel_indicator + '\n' + corrected_seq + '\n' + cigar_string +\
    #        '\n' + MD_seq + '\n' + aln.get_reference_sequence() + '\n' + aln.seq
    if (aln.is_reverse and aln.is_read1) or (not aln.is_reverse and aln.is_read2):
        ## always make seq from 5' to 3', position 0 == 5' of the fragment
        # end since RT comes from that side
        corrected_seq, indel_indicator = reverse_complement(corrected_seq[::-1]), indel_indicator[::-1]
    indel_indicator = np.array(list(indel_indicator))
    insertion_positions = np.where(indel_indicator=='I')[0]
    deletion_positions = np.where(indel_indicator=='D')[0]

    ## triplet pattern
    if insertion_positions.size != 0:
        indel_pattern_dict =  adding_indel(indel_pattern_dict,'Insertion',
                                           corrected_seq, insertion_positions)
    if deletion_positions.size != 0:
        indel_pattern_dict =  adding_indel(indel_pattern_dict,'Deletion',
                                           corrected_seq, deletion_positions)
    return indel_pattern_dict
    

def construct_df(indel_pattern_dict, filename):
    dfs = []
    for indel, pattern_dict in indel_pattern_dict.iteritems():
        df = pd.DataFrame({'pattern':pattern_dict.keys(),
                          'counts':pattern_dict.values()})\
            .assign(indel = indel)
        if len(df) > 0:
            dfs.append(df)
    if len(dfs)>0:
        return pd.concat(dfs, axis = 0) \
                .assign(samplename = filename)


    
def evaluate_aln(aln):
    pass_aln = not aln.is_supplementary and not aln.is_secondary
    pass_indel = 'I' in aln.cigarstring or 'D' in aln.cigarstring
    return pass_aln and pass_indel
    

def run_bam(bam_file):
    """ 
    Frame work for running bam file
    """
    filename = os.path.basename(bam_file).split('.')[0]
    print 'Running %s' %(filename)
    indel_pattern_dict = defaultdict(lambda: defaultdict(int))
    with pysam.Samfile(bam_file,'rb') as bam:
        for aln in bam:
            if not aln.is_unmapped and evaluate_aln(aln):
                indel_pattern_dict = analyze_alignment(aln, indel_pattern_dict)
    df = construct_df(indel_pattern_dict, filename)
    if df is not None:
        return df

            
    
def main():
    project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome'
    bam_path = project_path +  '/picard_results'
    indel_pattern_table_path = project_path + '/indel_pattern_tables'
    if not os.path.isdir(indel_pattern_table_path):
        os.makedirs(indel_pattern_table_path)
    tablename = indel_pattern_table_path + '/indel_pattern.tsv'
    figurename = tablename.replace('.tsv','.pdf')
    bam_files = glob.glob(bam_path + '/*.MarkDuplicate.bam')
    p = Pool(12)
    df = p.map(run_bam, bam_files)
    #df = map(run_bam, bam_files)
    p.close()
    p.join()
    df = filter(lambda x: x is not None,df)
    df = pd.concat(df, axis = 0)
    df.to_csv(tablename, index = False, sep='\t')
    print 'Written %s' %tablename
    #df = pd.read_table(tablename)
    #plot_indel(df, figurename)
    #print 'Written %s' %figurename
    

if __name__ == '__main__':
    main()
