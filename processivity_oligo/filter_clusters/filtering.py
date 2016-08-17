#!/usr/bin/env python

import pysam
from functools import partial
from multiprocessing import Pool
import os
import glob

def clusterCount(aln):
    name = aln.query_name
    return int(name.split('_')[-2])

def filterBam(result_path, min_read_cluster,bam_file):
    samplename = bam_file.split('/')[-1].split('.')[0]
    filtered_bam = result_path + '/' + samplename + '.filtered.bam'

    with pysam.Samfile(bam_file,'rb') as in_bam:
        with pysam.Samfile(filtered_bam,'wb',template = in_bam) as out_bam:
            [out_bam.write(aln) for aln in in_bam if clusterCount(aln) > min_read_cluster]

    print 'Written %s' %filtered_bam
    return 0

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/processivity'
    bam_path = project_path + '/bamFile'
    result_path = project_path + '/filtered_bams'
    min_read_cluster = 3

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
