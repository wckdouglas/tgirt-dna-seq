#!/usr/bin/env python

import pandas as pd
import numpy as np
from scipy.signal import savgol_filter
import os
import sys
import pyBigWig as pbw
import glob
from multiprocessing import Pool
from functools import partial
from itertools import izip
import argparse

def get_opt():
    chromosomes = range(1,23)
    chromosomes.extend(list('XY'))
    chromosomes = map(str, chromosomes)
    parser = argparse.ArgumentParser(description='Given a WPS file in bigwig format '+\
                                     'output peak coordinates in  bed file')
    parser.add_argument('-i', '--in_bigwig', help = 'Input bigWig', required=True)
    parser.add_argument('-o', '--out_bed', help = 'Output bed file name', required=True)
    parser.add_argument('-l','--length_type', help = 'short or long WPS?', default = 'Long',
                        choices = ['Short','Long'])
    parser.add_argument('-c','--chrom', help='Which chromosome?', choices = chromosomes, required=True)
    args = parser.parse_args()
    return args

def adjust_median(wps_array):
    """
    running median filter
    """
    rolling_median = pd.Series(wps_array)\
            .rolling(window = 1000)\
            .median()
    adjusted_wps = np.nan_to_num(wps_array - rolling_median)
    return adjusted_wps


def find_peak_region(wps):
    """
    given a wps array, output list of cooridinate of:
    start and end of peaks
    looking for the positions that across 0 in the de-noised signal
    """
    signs = np.sign(wps)
    signs[signs==0] = -1
    start = np.where(np.diff(signs)>0)[0]
    end = np.where(np.diff(signs)<0)[0]
    return start, end


def merge_peaks(peak_start, peak_end):
    new_start = []
    new_end = []
    tolerance_unprotected = 5
    i = 0
    while i < len(peak_start)-2:
        new_start.append(peak_start[i])
        j = i
        while j < len(peak_start)-2 and peak_start[j+1] - peak_end[j] <= tolerance_unprotected:
            j += 1
        new_end.append(peak_end[j])
        j += 1
        i = j
    new_start.append(peak_start[i])
    new_end.append(peak_end[i])
    return np.array(new_start), np.array(new_end)

def maxSubArray(ls):
    '''
    #https://gist.github.com/alabid/3734606
    '''
    if len(ls) == 0:
       raise Exception("Array empty") # should be non-empty

    runSum = maxSum = ls[0]
    i = 0
    start = finish = 0

    for j in range(1, len(ls)):
    	if ls[j] > (runSum + ls[j]):
            runSum = ls[j]
            i = j
        else:
            runSum += ls[j]

        if runSum > maxSum:
            maxSum = runSum
            start = i
            finish = j

    return start, finish


def pick_peak(above_median_starts, above_median_ends, sub_wps):
    '''
        from region that has 50 < size < 150,
        pick best peak (with maximum wps score)
    '''
    sub_wps = np.asarray(sub_wps)
    above_median_ends = np.asarray(above_median_ends)
    above_median_starts = np.asarray(above_median_starts)
    max_wps_array = np.array([sub_wps[s:e].max() for s, e in izip(above_median_starts, above_median_ends)])
    maximum_wps = np.where(max_wps_array == max_wps_array.max())[0]
    return above_median_starts[maximum_wps], above_median_ends[maximum_wps]


def adjust_peaks(wps, peak_start, peak_end):
    sub_wps = wps[peak_start:peak_end]
    median_sub_wps = np.median(sub_wps)
    adjusted_sub_wps = sub_wps - median_sub_wps
    return adjusted_sub_wps


def calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count, outFile, peak_size_filter = True):
    '''
        using peak start and end from wps array,
        find maximum sub array
        and export start and end from maximum subarray
        peak score determine from maximum wps score.
    '''
    adjusted_sub_wps = adjust_peaks(wps, peak_start, peak_end)
    above_median_starts, above_median_ends = find_peak_region(adjusted_sub_wps)

    if len(above_median_starts)>len(above_median_ends):
        above_median_ends = np.append(above_median_ends,len(adjusted_sub_wps))
    if not peak_size_filter:
        above_median_starts, above_median_ends =  pick_peak(above_median_starts, above_median_ends, adjusted_sub_wps)

    for above_median_start, above_median_end in izip(above_median_starts, above_median_ends):
        sub_peak_wps = adjusted_sub_wps[above_median_start:above_median_end]
        nucleosome_start , nucleosome_end = maxSubArray(sub_peak_wps)


        #adjust coordinate
        nucleosome_start, nucleosome_end = peak_start + above_median_start + np.array([nucleosome_start, nucleosome_end])
        nucleosome_center = int((nucleosome_start + nucleosome_end) /2)
        peak_center = (nucleosome_start + nucleosome_end)/2
        nucleosome_size = abs(nucleosome_end - nucleosome_start)
        if (peak_size_filter and 50 < nucleosome_size  < 150 ) or (not peak_size_filter and nucleosome_size > 5):
            peak_score = wps[nucleosome_start:nucleosome_end].max()
            peak_count += 1
            peak_name = '%s_peak%i' %(chromosome, peak_count)
            line = '\t'.join(map(str,[chromosome, nucleosome_start, nucleosome_end, peak_name, peak_score, '+', peak_center]))
            outFile.write(line+'\n')
    return peak_count




def write_short_peaks(wps, start, end, out_bed, chromosome):
    peak_count = 0
    for peak_start, peak_end in izip(start, end):
        peak_count += 1
        peak_name = '%s_peak%i' %(chromosome, peak_count)
        peak_score = wps[peak_start:peak_end].max()
        peak_center = long((peak_end - peak_start) /2)
        variables = map(str,[chromosome, peak_start, peak_end, peak_name, peak_score, '+', peak_center])
        line = '\t'.join(variables)
        out_bed.write(line+'\n')
    return peak_count


def write_long_peaks(wps, start, end, out_bed, chromosome):
    peak_count = 0
    start, end = merge_peaks(start, end)
    for peak_start, peak_end  in izip(start, end):
        peak_size = np.abs(peak_end - peak_start)
        if 40 <= peak_size <= 150:
            peak_count = calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count,
                                            out_bed, peak_size_filter = False)
        elif 150 < peak_size <= 450:
            peak_count = calling_long_peaks(chromosome, wps, peak_start, peak_end, peak_count,
                                            out_bed, peak_size_filter = True)
    return peak_count


def write_peaks(wps, peak_bed, length_type, chrom):
    with open(peak_bed, 'w') as out_bed:
        start, end = find_peak_region(wps)
        if length_type == 'Long':
            peak_count = write_long_peaks(wps, start, end, out_bed, chrom)
        elif length_type == 'Short':
            peak_count = write_short_peaks(wps, start, end, out_bed, chrom)
        else:
            sys.exit('Undetermined length type')
    print 'Written %i peaks to %s' %(peak_count, peak_bed)
    return 0


def process_bigwig(peak_bed, inputWig, chromosome, length_type):
    # get bigwig information
    filename = os.path.basename(inputWig)

    #print message
    print 'Running %s as %s for chrom: %s' %(filename, length_type, chromosome)

    # read in data
    bw = pbw.open(inputWig)
    chrom, length = bw.chroms().items()[0]
    assert chrom == chromosome, 'Wrong chromosomes'
    wps = np.array(bw.values(chrom,0,length))
    bw.close()
    wps = adjust_median(wps)
    if length_type == 'Long':
        wps = savgol_filter(wps, window_length = 21, polyorder=2)

    write_peaks(wps, peak_bed, length_type, chrom)
    return 0


def main():
    args = get_opt()
    bigwig_file = args.in_bigwig
    peak_bed = args.out_bed
    chromosome = args.chrom
    length_type = args.length_type
    process_bigwig(peak_bed, bigwig_file, chromosome, length_type)
    return 0


if __name__ == '__main__':
    main()
