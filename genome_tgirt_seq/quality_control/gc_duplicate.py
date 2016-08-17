#!/usr/bin/env python

import pandas as pd
from glob import glob
import os
from collections import defaultdict
from multiprocessing import Pool
from pybedtools import BedTool, set_tempdir
from pybedtools.cbedtools import Interval
from functools import partial
import numpy as np

def filterCount(frag):
    return int(frag.fields[-1])  > 1

def restructureFrag(frag):
    return Interval(chrom = frag.chrom,
            start = frag.start,
            end = frag.end,
            name = frag.fields[6],
            score = frag.fields[4],
            strand = frag.fields[3])

def runBed(ref_fa, result_path, bed_file):
    samplename = os.path.basename(bed_file).split('.')[0]
    result_file = result_path + '/' + samplename + '.tsv'
    print 'Running %s' %samplename
    df = BedTool(bed_file) \
        .groupby(g='1,2,3,6', o='count',c=1) \
        .filter(filterCount) \
        .nucleotide_content(fi = ref_fa) \
        .each(restructureFrag) \
        .moveto(result_file) \
        .to_dataframe(names = ['chrom','start','end','gc_per','count','strand'],
                dtype = {'chrom':str,'start':np.int32, 'end':np.int32,
                            'gc_per':np.float64, 'count':np.int32, 'strand':str})\
        .assign(length = lambda d: d['end'] - d['start']) \
        .pipe(lambda d: d[['gc_per','count','length']]) \
        .assign(gc_per = lambda d: np.array(d['gc_per']*100,dtype=np.int32))\
        .assign(samplename = samplename)
    df.to_csv(result_file, index=False,sep='\t')
    print 'Written %s' %result_file
    return df

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bed_path = project_path + '/bedFiles/sorted_bed'
    result_path = project_path + '/gc_length_duplicate'
    ref_path = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA'
    ref_fa = ref_path + '/genome_rDNA.fa'
    result_tsv = result_path + '/dup_gc_length.tsv'

    if not os.path.isdir(result_path):
        os.mkdir(result_path)
    set_tempdir(result_path)

    bed_files = glob(bed_path + '/*.bed')
    func = partial(runBed, ref_fa, result_path)
    p = Pool(12)
    dfs = p.map(func, bed_files)
    p.close()
    p.join()

    df = pd.concat(dfs,axis=0)
    df.to_csv(result_tsv,sep = '\t',index=False)
    print 'Written %s ' %result_tsv


if __name__ == '__main__':
    main()
