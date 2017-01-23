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
from wps_spacing import daniell_spectrum, r_spectrum
from collections import defaultdict


def extract_fft(wps_array):
    periodicity, intensity = daniell_spectrum(wps_array)
    #periodicity, intensity = r_spectrum(wps_array)
    usable_indices = (periodicity<=280) & (periodicity>=120)
    periodicity = periodicity[usable_indices]
    intensity = intensity[usable_indices]
    return periodicity, intensity


def find_tss_periodicity(wps):
    result = None
    wps = np.asarray(wps)
    signs = np.sign(wps)
    signs[signs==0] = -1
    peak_count = np.where(np.diff(signs)>0)[0]
    #if len(peak_count) > 2:
    result =  extract_fft(wps)
    return result


def make_tss_region(start, end, strand):
    if strand == '-':
        tss_start = end - 10000
        tss_end = end
    elif strand == '+':
        tss_start = start
        tss_end = start + 10000
    return tss_start, tss_end


def run_chrom_genes(chrom_genes, bw, gene_count, out):
    for count, gene in chrom_genes.iterrows():
        chrom = gene['chrom']
        start = gene['start']
        end = gene['end']
        strand = gene['strand']
        tss_start, tss_end = make_tss_region(start, end, strand)

        wps_array = bw.values(chrom, tss_start, tss_end)
        result = find_tss_periodicity(wps_array)
        if result is not None:
            line_info = '\t'.join(map(str,[gene['name'] ,gene['type'], gene['id']]))
            for period, intensity in zip(*result):
                out.write('{info}\t{period}\t{intensity}\n'.format(info = line_info,
                                                            period = period,
                                                            intensity = intensity))
        gene_count += 1
        if gene_count % 5000 == 0:
            print 'Parsed {gene_count} for {filename}'.format(gene_count = gene_count,
                                                        filename = out.name)
    return gene_count


def run_file(protein_df, out_path, bw_path, bw_prefix):
    filename = bw_prefix
    out_file = out_path + '/' + bw_prefix + '.bed'
    bw_prefix = bw_path + '/' + bw_prefix

    fail = 0
    gene_count = 0
    with open(out_file,'w') as out:
        out.write('name\ttype\tid\tperiodicity\tintensity\n')
        chromosomes = protein_df.chrom.unique()
        for chrom in chromosomes:
            print 'Running {filename} for chrom: {chrom}'.format(filename=filename, chrom=chrom)
            chrom_genes = protein_df.query("chrom == '%s'" %chrom)
            bw_file = '%s.%s.Long.bigWig' %(bw_prefix,chrom)
            bw = pbw.open(bw_file, 'r')
            gene_count = run_chrom_genes(chrom_genes, bw, gene_count, out)
            bw.close()
    return 0


def genes_to_mem(protein_bed):
    colnames = ['chrom','start','end','name','score','strand','type','id']
    gene_df = pd.read_table(protein_bed, names = colnames)
    print 'Read genes'
    return gene_df


def main():
    ref_path = os.environ['REF']
    protein_bed = ref_path + '/GRCh38/Bed_for_counts_only/protein.bed'
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path + '/bigWig_files'
    out_path = project_path + '/tss_periodicity'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    bw_files = glob.glob(bw_path + '/*merged*.bigWig')
    bw_prefix = set(map(lambda x: os.path.basename(x).split('.')[0], bw_files))
    chroms = range(1,23)
    chroms.extend(['X','Y'])
    chroms = map(str, chroms)
    protein_df = genes_to_mem(protein_bed)\
        .pipe(lambda d: d[np.in1d(d.chrom, chroms)])
    run_file_func = partial(run_file, protein_df, out_path, bw_path)
    p = Pool(6)
    p.map(run_file_func, list(bw_prefix))
    p.close()
    p.join()
    return 0


if __name__ == '__main__':
    main()
