#!/usr/bin/env python -u

import glob
import os
from multiprocessing import Pool
from functools import partial
import pyBigWig as pbw
import numpy as np
import pandas as pd
from genomeWPS import writeWig
import re
import sys


def read_bw(chrom, filename):
    bw = pbw.open(filename)
    chrom_len = bw.chroms()[chrom]
    values = np.array(bw.values(chrom, 0, chrom_len))
    return values


def merge_file(resultpath, df, args):
    chrom, length_type = args
    outfile = '%s/P1022_1113_13_mixed_unique.%s.%s.bigWig' %(resultpath, chrom, length_type)
    if os.path.isfile(outfile):
        sys.exit('Same name? %s' %outfile)
    filtered_df = df.query('chrom == "%s" & length_type == "%s"'\
                           %(chrom, length_type))
    filenames = np.array(filtered_df.filenames)
    readWPS = partial(read_bw, chrom)
    merged_wps = np.sum(map(readWPS, filenames), axis = 0)
    writeWig(merged_wps, outfile, chrom)
    print 'Merged: %s for chrom %s' %(length_type, chrom)
    return 0

def main():
    projectpath = os.environ['WORK'] + '/cdw2854/plasmaDNA'
    datapath = projectpath + '/genomeWPS/bigWig_files'
    resultpath = projectpath + '/genomeWPS/bigWig_files'


    # make dataframe
    filenames = glob.glob(datapath + '/P1*')
    files = map(os.path.basename,filenames)
    files = filter(lambda x: re.search('P1[013]',x), files)
    files = filter(lambda x: re.search('umi2id_unique',x), files)
    files = filter(lambda x: not re.search('mix',x), files)
    df = pd.DataFrame(pd.Series(files).str.split('.',4).tolist(),
                      columns=['prefix','chrom','length_type','filetype']) \
            .assign(filenames = map(lambda x: datapath+'/'+x, files))

    chromosomes = df.chrom.unique()
    length_type = df.length_type.unique()
    iterable = [(chrom,length) for chrom in chromosomes for length in length_type]

    # start merging
    mergeFile = partial(merge_file, resultpath, df)
    Pool(12).map(mergeFile, iterable)
    return 0

if __name__ == '__main__':
    main()
