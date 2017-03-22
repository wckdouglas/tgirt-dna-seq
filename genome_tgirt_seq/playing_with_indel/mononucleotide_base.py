
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
from multiprocessing import Pool
from itertools import izip

complementary = string.maketrans('ACTG','TGAC')

def keep_indel(map_position,cigar_seq, start, end):
    cigar_seq = np.array(list(cigar_seq))
    map_position = np.array(map_position)
    usable = (map_position >= start - 1) & (map_position <= end + 1)
    cigar = cigar_seq[usable]
    return ''.join(cigar)


def cigar_to_seq(cigar):
    """
    Input a cigar string (eg. 1S16M5I10M4S)
    output a line: SMMMMMMMMMMMMMMMMIIIIIMMMMMMMMMMSSSS
    """
    cigarNum = np.array(re.findall('[0-9]+',cigar),dtype='int64')
    cigarStr = np.array(re.findall('[A-Z]',cigar),dtype='string')
    usable = np.in1d(cigarStr,np.array(['S','M','I','D'],dtype='string'))
    cigarStr = cigarStr[usable]
    cigarNum = cigarNum[usable]
    cigarSeq = ''
    for s, n in zip(cigarStr, cigarNum):
        cigarSeq += int(n)*str(s)
    return cigarSeq

def calibrate_seq(cigar_seq, sequence, ref_positions):
    """
    making cigar seq and seq as same length
    with Deletions as '-'
    """
    new_sequence = ''
    new_pos = []
    new_cigar = ''

    acceptable_cigar = list('M')
    seq = iter(sequence)
    pos = iter(ref_positions)
    for cigar in cigar_seq:
        if cigar == 'S':
            seq.next()
        elif cigar == 'D':
            new_cigar += cigar
            new_pos.append(current_position + 1)
            new_sequence += '-'
        elif cigar == 'I':
            new_cigar += cigar
            new_pos.append(current_position)
            current_base = seq.next()
            new_sequence += current_base

        elif cigar == 'M':
            current_base = seq.next()
            current_position = pos.next()
            new_sequence += current_base
            new_pos.append(current_position)
            new_cigar += cigar
    return new_cigar, new_sequence, new_pos


def get_strand(aln):
    cDNA_reverse =  aln.is_read1 and aln.is_reverse or aln.is_read2 and not aln.is_reverse
    return '+' if cDNA_reverse else '-'

def running_mononuclotide_regions(line, bam, out):
    fields = line.strip().split('\t')
    seq_id, start, end, run_length, mononucleotide  = fields[0], long(fields[1]), long(fields[2]), fields[4], fields[-1]
    mono_indel_count = defaultdict(lambda: defaultdict(int))
    output_line = ''
    pos_aln_count, neg_aln_count = 0, 0
    for aln in bam.fetch(seq_id, start, end):
        strand = get_strand(aln)
        if strand == '-':
            neg_aln_count += 1
        else:
            pos_aln_count += 1
        if not aln.is_unmapped and re.search('I|D', aln.cigarstring):
            cigar_seq = cigar_to_seq(aln.cigarstring)
            cigar_seq, sequence, ref_positions = calibrate_seq(cigar_seq, aln.query_sequence, aln.get_reference_positions())
            cigar = keep_indel(ref_positions, cigar_seq, start, end)
            indel = re.findall('[ID]', cigar)
            indel = set(indel) # normalized read
            cigar_counter = Counter(indel)
            mono_indel_count[strand]['deletion'] += cigar_counter['D']
            mono_indel_count[strand]['insertion'] += cigar_counter['I']
    if neg_aln_count + pos_aln_count !=0:
        #print line, neg_aln_count, pos_aln_count
        #print mono_indel_count
        #print cigar_seq + '\n' + sequence + '\n', ref_positions
        #print start, end, cigar
        output_line = '{nucleotide}\t{run_len}\t'.format(nucleotide = mononucleotide, run_len = run_length) +\
                '{fwd_insert}\t{fwd_delete}\t'.format(
                            fwd_insert = mono_indel_count['+']['insertion'],
                    fwd_delete = mono_indel_count['+']['deletion']) +\
                '{rvs_insert}\t{rvs_delete}\t'.format(
                            rvs_insert = mono_indel_count['-']['insertion'],
                            rvs_delete = mono_indel_count['-']['deletion']) +\
                '{pos_count}\t{neg_count}'.format(pos_count = pos_aln_count, neg_count = neg_aln_count)
    return output_line


def make_index(bam_file):
    bam = pysam.Samfile(bam_file,'rb')
    try:
        bam.check_index()
        bam.close()
    except ValueError:
        bam.close()
        pysam.index(bam_file)
        print 'Indexed %s' %bam_file
    return 0

def analyze_bam_files(out_path, indel_table, bam_file):
    samplename = os.path.basename(bam_file).replace('.bam','')
    out_file = out_path + '/' + samplename + '.tsv'
    print 'Analyzing %s' %out_file
    header = 'mononucleotide\trun_length\tfwd_insertion\tfwd_deletion\t'+\
            'rev_insertion\trev_deletion\t'+\
            'positive_aln_count\tnegatve_aln_count'
    make_index(bam_file)
    with pysam.Samfile(bam_file, 'rb') as bam, \
            open(indel_table, 'r') as indel_file, \
            open(out_file,'w') as out:
        out.write(header + '\n')
        for line in indel_file:
            outline = running_mononuclotide_regions(line, bam, out)
            if outline != '':
                out.write(outline + '\n')
    print 'Finished %s' %out_file
    return 0


def main():
    ref_path = '/stor/work/Lambowitz/ref/Ecoli'
    indel_table = ref_path + '/k12_mg1655_high_indel.bed'
    project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome'
    bam_path = project_path + '/picard_results'
    out_path = project_path + '/base_indel_table'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    bam_files = glob.glob(bam_path + '/K12*UMI*clustered_family*.subsampled.bam')
    bam_files = filter(lambda x: re.search('[0-9].subsampled.bam',x), bam_files)

    analyzing_function = partial(analyze_bam_files, out_path, indel_table)
    p = Pool(24)
    p.map(analyzing_function, bam_files)
    #map(analyzing_function, bam_files)
    p.close()
    p.join()
    return 0


if __name__ == '__main__':
    main()
