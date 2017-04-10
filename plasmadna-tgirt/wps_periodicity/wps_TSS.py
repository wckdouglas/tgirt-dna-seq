#!/usr/bin/env python

import pyBigWig as pbw
import numpy as np
from scipy import fftpack
from itertools import izip
from functools import partial
import pandas as pd
from multiprocessing import Pool
import glob
import os
from wps_spacing import daniell_spectrum, r_spectrum, r_spec_pgram
from collections import defaultdict
from rpy2.rinterface import RRuntimeError
import sys
import re
import argparse


r_periodogram = r_spec_pgram()
def extract_fft(wps_array, args):
    if args.fft_type == 'R':
        periodicity, intensity = r_spectrum(wps_array, r_periodogram)
    else:
        periodicity, intensity = daniell_spectrum(wps_array)
    usable_indices = (periodicity<=280) & (periodicity>=120)
    periodicity = periodicity[usable_indices]
    intensity = intensity[usable_indices]
    return periodicity, intensity


def find_tss_periodicity(wps, args):
    result = None
    wps = np.asarray(wps)
    signs = np.sign(wps)
    signs[signs==0] = -1
    peak_count = np.where(np.diff(signs)>0)[0]
    #if len(peak_count) > 2:
    result =  extract_fft(wps, args)
    return result


def make_tss_region(start, end, strand):
    if strand == '-':
        tss_start = end - 10000
        tss_end = end
    elif strand == '+':
        tss_start = start
        tss_end = start + 10000
    return tss_start, tss_end


def make_whole_gene(start, end ,strand):
    return start, end


def run_gene(args, gene, out):
    chrom = gene['chrom']
    start = gene['start']
    end = gene['end']
    strand = gene['strand']
    if args.gene_regions == 'tss':
        tss_start, tss_end = make_tss_region(start, end, strand)
    else:
        tss_start, tss_end = make_whole_gene(start, end, strand)

    wps_array = bw.values(chrom, tss_start, tss_end)
    try:
        result = find_tss_periodicity(wps_array, args)
        if result is not None:
            line_info = '\t'.join(map(str,[gene['name'] ,gene['type'], gene['id']]))
            for period, intensity in zip(*result):
                out.write('{info}\t{period}\t{intensity}\n'.format(info = line_info,
                                                        period = period,
                                                        intensity = intensity))
    except RRuntimeError:
        print gene
        sys.exit()
    return 0


def run_file(args):
    '''
    For genes in a chromosome, get tss/gene body and run fourier transform on the gene
    '''
    bw = pbw.open(args.in_bigwig, 'r')
    gene_count = 0
    with open(args.out_bed,'w') as out:
        out.write('name\ttype\tid\tperiodicity\tintensity\n')
        for count, gene in chrom_genes.iterrows():
            run_gene(args, gene, out)
            if count % 1000 == 0:
                print 'Parsed {gene_count} for {filename} at {chrom}'.format(gene_count = gene_count,
                                                        filename = out.name,
                                                        chrom = args.chrom)
    bw.close()
    return out_file

def genes_to_mem(protein_bed):
    '''
    read in gene bed file
    '''
    colnames = ['chrom','start','end','name','score','strand','type','id']
    gene_df = pd.read_table(protein_bed, names = colnames) \
            .query('end-start > 1000')
    print 'Read genes'
    return gene_df

chroms = range(1,23)
chroms.extend(['X','Y'])
chroms = map(str, chroms)
def get_opt():
    parser = argparse.ArgumentParser(description='FFT WPS signal at genes')
    parser.add_argument('-i','--in_bigwig', required=True, help='BigWig file containing WPS signal')
    parser.add_argument('-c','--chrom',required=True, help = 'chromosome to analyze', choices = chroms)
    parser.add_argument('-o','--out_bed', required=True, help = 'Out bed file storing FFT intensity and freq at each row')
    parser.add_argument('-t','--fft_type', default = 'scipy', choices = ['scipy','R'], help='FFT algorithm to use (default: scipy)')
    parser.add_argument('-g','--gene_bed', required=True, help = 'Bed file storing genes')
    parser.add_argument('-r','--gene_regions', default='tss', choices = ['tss','gene_body'], help = 'gene regions to take for FFT')
    args = parser.parse_args()
    return args

def main():
    args = get_opt()
    protein_df = genes_to_mem(args.gene_bed)\
        .pipe(lambda d: d[np.in1d(d.chrom, chroms)])
    chrom_genes = protein_df.query("chrom == '%s'" %args.chrom)
    filename = os.path.basename(args.in_bigwig)
    print 'Running {filename} for chrom: {chrom}'.format(filename=filename, chrom=args.chrom)

    return 0


if __name__ == '__main__':
    main()
