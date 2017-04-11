#!/usr/bin/env python

import numpy as np
import fileinput
from pybedtools import BedTool
import os
import pandas as pd
import time
from multiprocessing import Pool
from functools import partial
from collections import defaultdict


def parseClosetLine(closestLine):
    fields = closestLine.fields
    center1 = int(fields[6])
    center2 = int(fields[13])
    distance = center2 - center1
    return distance

def closestPeak(bedpath, file1, file2, chromosome):
    print 'Running chromosome: %s' %chromosome
    bedFile1 = '%s/%s.%s.Long.bed' %(bedpath, file1, chromosome)
    bedFile2 = bedFile1.replace(file1,file2)
    if os.path.isfile(bedFile1) and os.path.isfile(bedFile2):
        bed1 = BedTool(bedFile1)
        bed2 = BedTool(bedFile2)
        distance_count = defaultdict(int)
        for line in BedTool(bed1).closest(bed2):
            center_distance = parseClosetLine(line)
            if -1000 < center_distance < 1000:
                distance_count[center_distance] += 1
        df = pd.DataFrame({'distance':distance_count.keys(), 
                           'count':distance_count.values()}) \
            .assign(chrom = chromosome)
        return df

def makeDir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
    return 0

def main():
    start = time.time()
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bedpath = projectpath + '/genomeWPS/bed_files'
    tablename = projectpath + '/figures/predictedNucleosomeDistance.tsv'
    file1 = 'P1022_1113_13_1016_mix_unique'
    file2 = 'SRR2130051_rmdup'
    chromosomes = map(str,np.arange(1,23))
    chromosomes = np.concatenate([chromosomes,['X','Y']])
    closestPeakFunc = partial(closestPeak, bedpath, file1, file2)
    p = Pool(24)
    dfs = p.map(closestPeakFunc, chromosomes)
    dfs = [df for df in dfs if df is not None]
    df = pd.concat(dfs)
    df.to_csv(tablename, index=False, sep = '\t')
    df = pd.read_table(tablename)
    print 'Writtten %s in %.3f min' %(tablename, np.true_divide(time.time() - start,60))
    return 0

if __name__ == '__main__':
    main()
