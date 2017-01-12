#!/usr/bin/env python

import pyBigWig as pbw
import glob
import os

bw_file_path = '/scratch/02727/cdw2854/plasmaDNA/genomeWPS/bigWig_files'
bw_files = glob.glob(bw_file_path + '/*bigWig')

for bw in bw_files:
    try:
        pbw.open(bw,'r')
    except:
        print '.'.join(os.path.basename(bw).split('.')[:2])
