#!/usr/bin/env python

import pyximport
import pysam
import numpy as np
from functools import partial
from multiprocessing import Pool
import os
import glob
import sys

def clusterCount(aln):
    name = aln.query_name
    cluster_count = int(name.split('_')[-2])
    return cluster_count

def filterBam(result_path, min_read_cluster, bam_file):

    samplename = bam_file.split('/')[-1].split('.')[0]
    filtered_bam = result_path + '/' + samplename + '.filtered.%i.bam' %(min_read_cluster)

    with pysam.Samfile(bam_file,'rb') as in_bam:
        with pysam.Samfile(filtered_bam,'wb',template = in_bam) as out_bam:
            [out_bam.write(aln) for aln in in_bam if clusterCount(aln) >= min_read_cluster]
    print 'Written %s' %filtered_bam
    os.system('samtools index %s' %filtered_bam)
    print 'Indexed %s' %filtered_bam
    return 0


def main():
    if len(sys.argv) != 2:
        sys.exit('Not enought argument\n'+
                '[usage] python %s <cluster member cutoff> \n')
    project_path = '/stor/work/Lambowitz/cdw2854/target-seq'
    bam_path = project_path + '/bamFiles'
    result_path = project_path + '/filtered_bams'
    min_read_cluster = int(sys.argv[1])


    if not os.path.isdir(result_path):
        os.mkdir(result_path)

    bam_files = glob.glob(bam_path + '/*bam')

    func = partial(filterBam,result_path, min_read_cluster)

    p = Pool(12)
    p.map(func, bam_files)
    p.close()
    p.join()

if __name__ == '__main__':
    main()
