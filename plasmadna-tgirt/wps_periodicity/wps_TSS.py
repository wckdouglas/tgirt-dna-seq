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
from wps_spacing import highest_periodicity


def find_tss_periodicity(bw_prefix, chrom, tss_start, tss_end):
    result = (0,0)
    bw_file = '%s.%s.Long.bigWig' %(bw_prefix,chrom)
    try:
        bw = pbw.open(bw_file, 'r')
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


chroms = np.array(np.arange(1,23),dtype='string')
def run_file(protein_bed, out_path, bw_path, bw_prefix):
    out_file = out_path + '/' + bw_prefix + '.bed'
    bw_prefix = bw_path + '/' + bw_prefix 

    preiodicity_func = partial(find_tss_periodicity, bw_prefix)
    with open(out_file,'w') as out:
        out.write('chrom\tstart\tend\tname\tscore\tstrand\ttype\tid\tperiodicity\tintensity\n')
        for bed_record in open(protein_bed, 'r'):
            fields = bed_record.split('\t')
            chrom = fields[0]
            start = long(fields[1])
            if chrom in chroms:
                end = long(fields[2])
                strand = fields[5]
            
                if strand == '-':
                    tss_start, tss_end  = end - 2000, end + 2000
                elif strand == '+':
                    tss_start, tss_end = start - 2000, start + 2000

                max_periodicity, max_intensity = preiodicity_func(chrom, tss_start, tss_end)
                out.write('%s\t%s\t%s\n' %(bed_record.rstrip(), str(max_periodicity),str(max_intensity)))
    print 'written: %s' %out_file
    return 0

def main():
    ref_path = os.environ['REF']
    protein_bed = ref_path + '/GRCh38/Bed_for_counts_only/protein.bed'
    project_path =  '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
    bw_path = project_path + '/bigWig_files'
    out_path = project_path + '/tss_periodicity'
    if not os.path.isdir(out_path):
        os.mkdir(out_path)
    bw_prefix = ['PD_merged', 'SRR2130052']
    run_file_func = partial(run_file, protein_bed, out_path, bw_path)
    Pool(2).map(run_file_func, bw_prefix)


if __name__ == '__main__':
    main()
