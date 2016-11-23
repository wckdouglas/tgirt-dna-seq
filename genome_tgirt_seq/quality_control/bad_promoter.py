#!/usr/bin/env python

from pybedtools import BedTool
import pandas as pd
import glob
from multiprocessing import Pool
import os
from functools import partial

def get_coverage(promoter_bed, genome, window, bam_file):
    name = bam_file.split('/')[-1].split('.')[0]
    print 'Running %s' %name
    df = BedTool(promoter_bed) \
        .slop(b=window, g = genome)\
        .coverage(d = True, b = bam_file)\
        .to_dataframe(names = ['chrom','start','end','name',
                                'score','strand','ratio','position','coverage']) \
        .assign(samplename = name) \
        .drop(['chrom','start','end','score','strand'],axis=1)
    print 'Finished making table: %s' %name
    return df


def main():
    window = 1900
    work_path = os.environ['WORK']
    ref_path = os.environ['REF']
    bam_path = work_path + '/cdw2854/genomeDNA/clustered_map'
    hg19 = ref_path + '/hg19'
    bad_promoter = hg19 + '/annotations/reduced_bad_promoters.bed'
    genome = hg19 + '/Sequence/WholeGenomeFasta/genome.txt'
    out_table_name = bam_path + '/bad_promoter_coverage.tsv'
    bam_files = glob.glob(bam_path + '/*.bam')
    cov_func = partial(get_coverage, bad_promoter, genome, window)
    p = Pool(12)
    dfs = p.map(cov_func, bam_files)
    p.close()
    p.join()
    pd.concat(dfs, axis=0) \
        .to_csv(out_table_name, sep='\t', index=False)
    print 'Written %s' %out_table_name

if __name__ == '__main__':
    main()
