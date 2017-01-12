
#!/usr/bin/env python

import pyBigWig as pbw
import numpy as np
from itertools import izip
from functools import partial
import pandas as pd
from multiprocessing import Pool
import glob
import os
import sys
from spectrum import pdaniell
from statsmodels.tsa.filters.filtertools import recursive_filter
from scipy.signal import detrend


def shift_array(signal):
    signal = np.array(signal)
    signal = np.append(signal[:300], signal)
    return signal

filter_freq = 1.0/np.arange(5,100,4)
def recursive_filter_function(signal):
    signal = shift_array(signal)
    filter_length = len(filter_freq)
    signal = recursive_filter(signal, filter_freq)[300:]
    return signal

def demean(arr):
    return arr-arr.mean()

def daniell_spectrum(signal, sample_rate):
    filtered_signal = detrend(signal, type='linear')
    filtered_signal = recursive_filter_function(filtered_signal)
    filtered_signal = demean(filtered_signal)
    # detrending and smoothing before spectrogram
    p = pdaniell(filtered_signal, 2, detrend='linear', sampling = sample_rate)
    p()
    psd = p.psd
    freq = p.frequencies()
    return 1/np.array(freq), np.array(psd)

def highest_periodicity(wps_array):
    periodicity, intensity = daniell_spectrum(wps_array)
    usable_indices = (periodicity<=280) & (periodicity>=120)
    periodicity = periodicity[usable_indices]
    intensity = intensity[usable_indices]
    argmax = np.argmax(intensity)
    max_periodicity = periodicity[argmax]
    max_intensity = intensity[argmax]
    return max_periodicity, max_intensity

def max_period_from_chrom(bw, chromosome, input_arg):
    '''
    From a bigWig value array
    using start and end to extract desired region,
    run FFT to get frequency domain and corersponding intensity
    '''
    bin_count, (start, end) = input_arg
    wps_array = bw.values(chromosome, start, end)

    signs = np.sign(wps_array)
    signs[signs==0] = -1
    peak_count = np.where(np.diff(signs)>0)[0]
    max_periodicity, max_intensity = 0, 0
    if len(peak_count) > 20:
        max_periodicity, max_intensity = highest_periodicity(wps_array)

    if bin_count % 10000 == 0:
        print 'Analyzed %i bin' %bin_count
    bin_name = '%s_%i' %(chromosome, bin_count)
    return (chromosome, start, end, bin_name , max_periodicity, max_intensity)

def analyze_file(samplename, bw_name):
    samplename = os.path.basename(bw_name).split('.')[0]
    chromosome = bw_name.split('.')[1]
    print 'Reading chromosome: %s from %s' %(chromosome, samplename)
    bw = pbw.open(bw_name,'r')
    chromosome_info = bw.chroms()
    window_size = 5000
    chrom_length = chromosome_info[chromosome]
    no_of_bins = chrom_length / window_size
    start_positions = np.linspace(0, chrom_length, no_of_bins)
    start_positions = np.array(start_positions, dtype=np.int64)
    end_positions = np.roll(start_positions,-1)[:-1]
    get_periodicity = partial(max_period_from_chrom, bw, chromosome)

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
    '''
    This program used chromosome splitted bigWig files to calculate nucleosome spacing.
    *** Long fragments WPS
    '''
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path  + '/bigWig_files'
    periodicity_path = project_path + '/periodicity_tables'
    if not os.path.isdir(periodicity_path):
        os.mkdir(periodicity_path)
    bw_files = glob.glob(bw_path + '/*.Long.bigWig')
    samplenames = set(map(lambda x: os.path.basename(x).split('.')[0], bw_files))
    read_sample = partial(analyze_sample, bw_path, periodicity_path)
    map(read_sample, iter(samplenames))
    return 0


if __name__ == '__main__':
    main()
