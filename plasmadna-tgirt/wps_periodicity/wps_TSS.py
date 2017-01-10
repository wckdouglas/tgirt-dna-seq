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
from wps_spacing import daniell_spectrum


def extract_fft(wps_array):
    periodicity, intensity = daniell_spectrum(wps_array, 0.1)
    usable_indices = (periodicity<=280) & (periodicity>=120)
    periodicity = periodicity[usable_indices]
    intensity = intensity[usable_indices]
    return periodicity, intensity


def find_tss_periodicity(bw_prefix, chrom, tss_start, tss_end):
    result = None
    bw_file = '%s.%s.Long.bigWig' %(bw_prefix,chrom)
    bw = pbw.open(bw_file, 'r')
    wps_array = bw.values(chrom, tss_start, tss_end)
    signs = np.sign(wps_array)
    signs[signs==0] = -1
    peak_count = np.where(np.diff(signs)>0)[0]
    if len(peak_count) > 2:
        result =  extract_fft(wps_array)
    bw.close()
    return result


def make_tss_region(start, end, strand):
    if strand == '-':
        tss_start = end - 10000
        tss_end = end
    elif strand == '+':
        tss_start = start
        tss_end = start + 10000
    return tss_start, tss_end

def run_file(protein_bed, out_path, bw_path, chroms, bw_prefix):
    filename = bw_prefix
    print 'Running {filename}'.format(filename=filename)
    out_file = out_path + '/' + bw_prefix + '.bed'
    bw_prefix = bw_path + '/' + bw_prefix

    preiodicity_func = partial(find_tss_periodicity, bw_prefix)
    fail = 0
    with open(out_file,'w') as out, open(protein_bed, 'r') as genes:
        out.write('name\ttype\tid\tperiodicity\tintensity\n')
        for gene_count, bed_record in enumerate(genes):
            fields = bed_record.rstrip().split('\t')
            chrom = fields[0]
            if chrom in chroms:
                start = long(fields[1])
                end = long(fields[2])
                strand = fields[5]
                tss_start, tss_end = make_tss_region(start, end, strand)

                result = preiodicity_func(chrom, tss_start, tss_end)
                if result is not None:
                    line_info = fields[3] +'\t' + fields[6] + '\t'+fields[7]
                    for period, intensity in zip(*result):
                        out.write('{info}\t{period}\t{intensity}\n'.format(info = line_info,
                                                                      period = period,
                                                                      intensity = intensity))
            if gene_count % 5000 == 0:
                print 'Parsed {gene_count} for {filename}'.format(gene_count = gene_count,
                                                                   filename = filename)
    return 0


def main():
    ref_path = os.environ['REF']
    protein_bed = ref_path + '/GRCh38/Bed_for_counts_only/protein.bed'
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path + '/bigWig_files'
    out_path = project_path + '/tss_periodicity'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    bw_files = glob.glob(bw_path + '/*.bigWig')
    bw_prefix = set(map(lambda x: os.path.basename(x).split('.')[0], bw_files))
    chroms = range(1,23)
    chroms.extend(['X','Y'])
    chroms = map(str, chroms)
    run_file_func = partial(run_file, protein_bed, out_path, bw_path, chroms)
    p = Pool(24)
    p.map(run_file_func, list(bw_prefix))
    p.close()
    p.join()
    return 0


if __name__ == '__main__':
    main()
