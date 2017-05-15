
#!/usr/bin/env python

import numpy as np
from collections import Counter, defaultdict
from pyfaidx import Fasta
import pysam
import glob
import re
from functools import partial
import os
import string
from itertools import izip
import argparse
from operator import itemgetter
from mononucleotide_base import make_index, get_strand, cigar_to_seq
from multiprocessing import Pool
from dmisc.sam_action import MDToSeq
from numba import jit

complementary = string.maketrans('ACTG','TGAC')
def reverse_complement(kmer):
    return kmer.translate(complementary)[::-1]


def keep_indel(map_position,cigar_seq, seq, start, end):
    cigar_seq = np.array(list(cigar_seq))
    map_position = np.array(map_position)
    usable = (map_position >= start-1) & (map_position <= end +1)
    cigar = cigar_seq[usable]
    seq = np.array(list(seq))
    seq = seq[usable]
    return ''.join(cigar), ''.join(seq)

@jit()
def extract_indel(sub_cigar, sub_seq, mononucleotide):
    deletions = ''
    insertions = ''
    for c, b in zip(sub_cigar, sub_seq):
        if c == 'D':
            deletions += mononucleotide
        elif c == 'I':
            insertions += b
    return insertions, deletions


@jit()
def calibrate_seq(cigar_seq, sequence, md_seq, ref_positions):
    """
    making cigar seq and seq as same length
    with Deletions as '-'
    """
    new_sequence = ''
    new_pos = []
    new_cigar = ''
    new_md = ''

    seq = iter(sequence)
    pos = iter(ref_positions)
    md = iter(md_seq)
    current_position = 0
    for cigar in cigar_seq:
        if cigar == 'S':
            seq.next()
        elif cigar == 'D':
            new_cigar += cigar
            new_pos.append(current_position + 1)
            new_sequence += '-'
            new_md += md.next()

        elif cigar == 'I':
            new_cigar += cigar
            new_pos.append(current_position)
            current_base = seq.next()
            new_sequence += current_base
            new_md += '+'

        elif cigar == 'M':
            current_base = seq.next()
            current_position = pos.next()
            new_sequence += current_base
            new_pos.append(current_position)
            new_cigar += cigar
            new_md += md.next()
    return new_cigar, new_sequence, new_pos, new_md


def running_mononuclotide_regions(inline, bam, out):
    '''
    For each bedline (representing each homopolymer runs region)
    count how many reads mapped to fwd strand and rvs strand
    for each alignment with insertion or deletions, extract the indel bases
    '''
    fields = inline.strip().split('\t')
    seq_id, start, end, run_length, mononucleotide  = fields[0], long(fields[1]), long(fields[2]), fields[4], fields[-1]
    mono_indel_count = defaultdict(lambda: defaultdict(list))
    output_line = ''
    pos_aln_count, neg_aln_count = 0, 0
    strands = ['+','-']
    indel_types = ['deletions','insertions']
    for aln in bam.fetch(seq_id, start, end):
        strand = get_strand(aln)
        if strand == '-':
            neg_aln_count += 1
        else:
            pos_aln_count += 1
        if not aln.is_unmapped and re.search('I|D', aln.cigarstring):
            cigar_seq = cigar_to_seq(aln.cigarstring)
            md_seq = MDToSeq(aln.get_tag('MD'))
            cigar_seq, sequence, ref_positions, md_seq = calibrate_seq(cigar_seq, aln.query_sequence, md_seq, aln.get_reference_positions())
            sub_cigar, sub_seq = keep_indel(ref_positions, cigar_seq, sequence, start, end)
            insertions, deletions = extract_indel(sub_cigar, sub_seq, mononucleotide)
            if strand == '-':
                deletions = reverse_complement(deletions)
                insertions = reverse_complement(insertions)
                mono_indel_count[strand]['deletions'].append(deletions)
                mono_indel_count[strand]['insertions'].append(insertions)
            elif strand == '+':
                mono_indel_count[strand]['deletions'].append(deletions)
                mono_indel_count[strand]['insertions'].append(insertions)
    outlines = ''
    for strand in strands:
        coverage = pos_aln_count if strand == '+' else neg_aln_count
        ref_base = mononucleotide if strand == '+' else mononucleotide.translate(complementary)
        for indel_type in indel_types:
            indel_pattern = ','.join(mono_indel_count[strand][indel_type])
            indel_pattern = indel_pattern.strip(',')
            if coverage > 0:
                coverage = str(coverage)
                outline = '\t'.join([run_length, ref_base,strand, coverage, indel_type, indel_pattern])
                outlines +=  outline + '\n'
    return outlines

def analyze_bam_files(out_path, indel_table, bam_file):
    samplename = os.path.basename(bam_file).replace('.bam','')
    out_file = out_path + '/' + samplename + '.tsv'
    print 'Analyzing %s' %out_file
    header = 'run_length\tmononucleotide\tstrand\tcoverage\tindel\tpatterns'
    make_index(bam_file)
    with pysam.Samfile(bam_file, 'rb') as bam, \
            open(indel_table, 'r') as indel_file, \
            open(out_file,'w') as out:
        out.write(header + '\n')
        for line in indel_file:
            outline = running_mononuclotide_regions(line, bam, out)
            if outline != '':
                out.write(outline)
    print 'Finished %s' %out_file
    return 0


def main():
    ref_path = '/stor/work/Lambowitz/ref/Ecoli'
    indel_table = ref_path + '/k12_mg1655_high_indel.bed'
    project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome'
    bam_path = project_path + '/picard_results'
    out_path = project_path + '/base_insertion_table'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    bam_files = glob.glob(bam_path + '/K12*UMI*clustered_family*.subsampled.bam')
    bam_files = filter(lambda x: re.search('[0-9].30X.subsampled.bam',x), bam_files)

    analyzing_function = partial(analyze_bam_files, out_path, indel_table)
    p = Pool(24)
    p.map(analyzing_function, bam_files)
    #map(analyzing_function, bam_files)
    p.close()
    p.join()
    return 0




if __name__ == '__main__':
    main()
