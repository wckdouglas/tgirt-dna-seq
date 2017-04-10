#!/usr/bin/env pyrhon

import glob
import sys
import numpy as np
import os
from itertools import product

project_path = os.environ['WORK'] + '/cdw2854/plasmaDNA'
project_path = os.environ['SCRATCH'] + '/plasmaDNA'
in_path = project_path+ '/bedFiles/chrom_split_bed/demultiplexed'
out_path = in_path + '/merged'
if not os.path.isdir(out_path):
    os.mkdir(out_path)
bed_files = glob.glob(in_path + '/*.bed')
bed_files = np.array(bed_files)
filenames = map(os.path.basename, bed_files)
samplename = map(lambda x: x.split('.')[0], filenames)
samplename = np.array(samplename)

for sn in np.unique(samplename):
    correct_files = (samplename == sn)
    get_files = bed_files[correct_files]
    files = ' '.join(get_files)
    command = 'cat %s > %s/%s.bed' %(files, out_path, sn)
    print command
