#!/bin/env python

import glob
import subprocess
import re
import os

datapath='/scratch/02727/cdw2854/plasmaDNA/splittedFastq'
fastqfiles = glob.glob(datapath+ '/*gz')
for f in fastqfiles:
    if not f.endswith('R1_001.fastq.gz') and not f.endswith('R2_001.fastq.gz'):
        basename = os.path.basename(f)
        name = basename.split('.')[0]
        suffix = '.'.join(basename.split('.')[1:])
        splittedNum = name.split('_')[-1]
        samplename = name.split('_')[0]
        readInfo = '_'.join(name.split('_')[1:-1])
        newname = '_'.join([samplename,splittedNum,readInfo]) + '.' + suffix 
        rename = 'mv %s %s/%s' %(f,datapath,newname)
        print rename
        subprocess.call(rename,shell=True)
