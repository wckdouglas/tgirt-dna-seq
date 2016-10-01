#!/usr/bin/env python

import pyximport
import pysam
import numpy as np
pyximport.install(setup_args={
        'extra_link_args':pysam.get_libraries(),
        'include_dirs':pysam.get_include(),
        'define_macros':pysam.get_defines()
        })
from functools import partial
from multiprocessing import Pool
import os
import glob
from cluster_filter import filterBam
import sys

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
