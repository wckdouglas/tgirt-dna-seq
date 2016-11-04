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


def highest_periodicity(bw, chromosome, input_arg):
    bin_count, (start, end) = input_arg
    wps_array = bw.values(chromosome, start, end)

    signs = np.sign(wps_array)
    signs[signs==0] = -1
    peak_count = np.where(np.diff(signs)>0)[0]
    max_periodicity, max_intensity = 0, 0
    if len(peak_count) > 50:
        periodicity, intensity = fft(wps_array)
        usable_indices = periodicity<500
        periodicity = periodicity[usable_indices]
        intensity = intensity[usable_indices]

        argmax = np.argmax(intensity)
        max_periodicity = periodicity[argmax]
        max_intensity = intensity[argmax]

    if bin_count == 1000:
        print 'Analyzing %i bin' %bin_count
    bin_name = '%s_%i' %(chromosome, bin_count)
    return (chromosome, start, end, bin_name , max_periodicity, max_intensity)

def analyze_file(samplename, bw_name):
    samplename = os.path.basename(bw_name).split('.')[0]
    bw = pbw.open(bw_name,'r')
    chromosome = bw_name.split('.')[1]
    print 'Reading %s from %s' %(chromosome, samplename)
    chromosome_info = bw.chroms()
    window_size = 10000
    chrom_length = chromosome_info[chromosome]
    no_of_bins = chrom_length / window_size
    start_positions = np.linspace(0, chrom_length, no_of_bins)
    start_positions = np.array(start_positions, dtype=np.int64)
    end_positions = np.roll(start_positions,-1)[:-1]
    get_periodicity = partial(highest_periodicity, bw, chromosome)

    iterator =  enumerate(izip(start_positions, end_positions))
    results = map(get_periodicity, iterator)
    df = pd.DataFrame(results, columns = ['chrom','start','end','name',
                                          'periodicity', 'intensity']) \
        .assign(samplename = samplename)
    bw.close()
    return df

def analyze_sample(bw_path, periodicity_path, samplename):
    print 'reading %s' %samplename
    file_regex = bw_path + '/' + samplename + '*Long.bigWig'
    bw_files = glob.glob(file_regex)
    bw_files = filter(lambda x: os.path.getsize(x) > 0, bw_files)
    read_file = partial(analyze_file, samplename)
    p = Pool(24)
    dfs = p.map(read_file, bw_files)
    p.close()
    p.join()
    df = pd.concat(dfs, axis = 0)
    periodicity_table = periodicity_path + '/' + samplename + '.tsv'
    df.to_csv(periodicity_table, index=False, sep ='\t')
    print 'Written %s' %periodicity_path
    return 0

def main():
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path + '/bigwig_files'
    periodicity_path = project_path + '/periodicity_tables'
    print periodicity_path
    if not os.path.isdir(periodicity_path):
        os.mkdir(periodicity_path)
    bw_files = glob.glob(bw_path + '/*.Long.bigWig')
    samplenames = set(map(lambda x: os.path.basename(x).split('.')[0], bw_files))
    read_sample = partial(analyze_sample, bw_path, periodicity_path)
    map(read_sample, iter(samplenames))
    return 0


if __name__ == '__main__':
    main()
