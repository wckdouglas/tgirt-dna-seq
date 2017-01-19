#!/usr/bin/env python -u

import glob
import os
from multiprocessing import Pool
from functools import partial
import pyBigWig as pbw
import numpy as np
import pandas as pd
from genomeWPS import writeWig


def read_bw(chrom, filename):
    bw = pbw.open(filename)
    chrom_len = bw.chroms()[chrom]
    values = np.array(bw.values(chrom, 0, chrom_len))
    return values


def merge_file(resultpath, df, args):
    chrom, prefix, length_type, process = args
    outfile = '%s/%s_%s_merged.%s.%s.bigWig' %(resultpath, prefix, process, chrom, length_type)
    if os.path.isfile(outfile):
        sys.exit('Same name? %s' %outfile)
    filtered_df = df.query('chrom == "%s" & length_type == "%s" & prefix == "%s" & process=="%s"'\
                           %(chrom, length_type, prefix,process))
    filenames = np.array(filtered_df.filename)
    readWPS = partial(read_bw, chrom)
    merged_wps = np.sum(map(readWPS, filenames), axis = 0)
    writeWig(merged_wps, outfile, chrom)
    print 'Merged: %s, %s, %s for chrom %s' %(prefix, length_type, process, chrom)
    return 0

def assign_prep(x):
    if 'umi2id_unique' in x:
        return 'UMI'
    elif 'clustered_rmdup' in x:
        return 'clustered_rmdup'
    elif 'clustered' in x:
        return 'clustered'
    elif 'rmdup' in x:
        return 'rmdup'
    else:
        return 'no_processing'

def assign_prefix(x):
    return x.split('_')[0].split('-')[0]


def main():
    projectpath = os.environ['WORK'] + '/cdw2854/plasmaDNA'
    datapath = projectpath + '/genomeWPS/bigWig_files'
    resultpath = projectpath + '/genomeWPS/bigWig_files'


    # make dataframe
    filenames = glob.glob(datapath + '/P1*')
    files = pd.Series(map(os.path.basename,filenames))
    df = pd.DataFrame(files.str.split('.',4).tolist(),columns=['prefix','chrom','length_type','filetype'])\
        .assign(filename = filenames) \
        .assign(process = lambda d: map(assign_prep, d.prefix))\
        .assign(prefix = lambda d: map(assign_prefix, d.prefix))

    chromosomes = df.chrom.unique()
    prefix = df.prefix.unique()
    length_type = df.length_type.unique()
    processes = df.process.unique()
    iterable = [(chrom,sample,length, process) for chrom in chromosomes for sample in prefix for length in length_type for process in processes]

    # start merging
    mergeFile = partial(merge_file, resultpath, df)
    Pool(12).map(mergeFile, iterable)
    return 0

if __name__ == '__main__':
    main()
