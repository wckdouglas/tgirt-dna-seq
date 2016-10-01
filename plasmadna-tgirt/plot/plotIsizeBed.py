#!/usr/bin/env python -u

import numpy as np
import sys
import time
import glob
import pandas as pd
import os
from multiprocessing import Pool
import plotIsizeBam

def getIsize(bed_file, samplename):
    print 'Analyzing %s ' %samplename
    with open(bed_file,'r') as bed:
        isize_array = [line.split('\t')[-2] for line in bed]
    print 'Extacted %i alignments from %s' %(len(isize_array), samplename)
    return np.array(isize_array,dtype=np.int32)

def main(bed_file):
    samplename = os.path.basename(bed_file).split('.')[0].split('_')[0]
    print 'Running %s ' %samplename
    isize_array = getIsize(bed_file, samplename)
    df = pd.DataFrame({'isize':isize_array[isize_array < 600]})\
	.assign(count = 1) \
	.groupby(['isize']) \
	.agg({'count':np.sum})\
	.reset_index()  \
	.assign(samplename = samplename) 
    print 'Finished %s ' %samplename
    return df

if __name__ == '__main__':
    start = time.time()
    project_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    datapath = project_path + '/bedFiles'
    figurepath = project_path + '/figures'
    figurename = figurepath + '/isizeBedTable.png'
    tablename = datapath + '/isizeBedTable.tsv'
    bed_files = glob.glob(datapath + '/*bed')
    p = Pool(24)
    dfs = p.map(main, bed_files)
    p.close()
    p.join()
    df = pd.concat(dfs)
    df.to_csv(tablename, sep='\t', index_label = False, index=False)
    df = pd.read_table(tablename,sep='\t')\
	.groupby(['samplename']) \
	.apply(plotIsizeBam.normCount) \
	.reset_index() 
    plotIsizeBam.plotFigure(figurename, df)
    print 'Time lapsed: %.3f sec' %(time.time() - start)
