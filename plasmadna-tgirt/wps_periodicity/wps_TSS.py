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


def fft(array):
    sample_size = len(array)
    half_size = sample_size/2
    intensity = fftpack.fft(array)
    intensity = abs(intensity)**2
    freq = fftpack.fftfreq(sample_size)
    periodicity = 1/(freq)
    return periodicity[:half_size], intensity[:half_size]


def highest_periodicity(wps_array):
    periodicity, intensity = fft(wps_array)
    usable_indices = periodicity<500
    periodicity = periodicity[usable_indices]
    intensity = intensity[usable_indices]

    argmax = np.argmax(intensity)
    max_periodicity = periodicity[argmax]
    max_intensity = intensity[argmax]
    return max_periodicity, max_intensity


def find_tss_periodicity(bw_prefix, chrom, tss_start, tss_end):
    result = (0,0)
    try:
        bw = pbw.open('%s.%s.Long.bigWig' %(bw_prefix,chrom))
        wps_array = bw.values(chrom, tss_start, tss_end)
        signs = np.sign(wps_array)
        signs[signs==0] = -1
        peak_count = np.where(np.diff(signs)>0)[0] 
        if len(peak_count) > 15:
            result =  highest_periodicity(wps_array)
        bw.close()
    except RuntimeError:
        pass
    return result


def main():
    ref_path = os.environ['REF']
    protein_bed = ref_path + '/GRCh38/Bed_for_counts_only/protein.bed'
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path + '/bigwig_files'
    bw_prefix = bw_path + '/PD_merged'
    bw_prefix = bw_path + '/SRR2130051'

    preiodicity_func = partial(find_tss_periodicity, bw_prefix)
    for bed_record in open(protein_bed, 'r'):
        fields = bed_record.split('\t')
        chrom = fields[0]
        start = long(fields[1])
        end = long(fields[2])
        strand = fields[5]
        interval =2
        
        if strand == '-':
            tss_start, tss_end  = end - 2000, end + 2000
        elif strand == '+':
            tss_start, tss_end = start - 2000, start + 2000

        max_periodicity, max_intensity = preiodicity_func(chrom, tss_start, tss_end)
        print '%s\t%s\t%s' %(bed_record.rstrip(), str(max_periodicity),str(max_intensity))
    return 0


if __name__ == '__main__':
    main()
