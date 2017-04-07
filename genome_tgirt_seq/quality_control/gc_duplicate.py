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
            name = frag.fields[3],
            score = frag.fields[7],
            strand = frag.fields[5])

def runBed(ref_fa, result_path, bed_file):
    samplename = os.path.basename(bed_file).split('.')[0]
    result_file = result_path + '/' + samplename + '.tsv'
    regular_chrom = map(str,range(1,20))
    print 'Running %s' %samplename
    df = BedTool(bed_file) \
        .filter(lambda f: f.chrom in regular_chrom )\
        .nucleotide_content(fi = ref_fa) \
        .each(restructureFrag) \
        .moveto(result_file) \
        .to_dataframe(names = ['chrom','start','end','cluster_id','gc_pct','strand'],
                dtype = {'chrom':str,'start':np.int32, 'end':np.int32,
                            'gc_pct':np.float64, 'cluster_id':str, 'strand':str})\
        .assign(length = lambda d: d['end'] - d['start']) \
        .assign(count = lambda d: np.array(d.cluster_id.str.split('_').str.get(1),dtype=np.int32)) \
        .pipe(lambda d: d[['gc_pct','count','length']]) \
        .assign(gc_per = lambda d: np.array(d['gc_pct']*100,dtype=np.int32))\
        .assign(samplename = samplename)
    df.to_csv(result_file, index=False,sep='\t')
    print 'Written %s' %result_file
    return df

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bed_path = project_path + '/bedFiles'
    result_path = project_path + '/gc_length_duplicate'
    ref_path = '/stor/work/Lambowitz/ref/GRCh38/hg38_rDNA'
    ref_fa = ref_path + '/genome_rDNA.fa'
    result_tsv = result_path + '/dup_gc_length.tsv'

    if not os.path.isdir(result_path):
        os.mkdir(result_path)
    set_tempdir(result_path)

    bed_files = glob(bed_path + '/*_unique.bed')
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
