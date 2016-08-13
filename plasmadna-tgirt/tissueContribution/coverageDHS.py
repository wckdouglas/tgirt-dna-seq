#!/usr/bin/env python

from pybedtools import BedTool
import os
import glob
from multiprocessing import Pool
import re

def coverage(args):
    bedFile, resultpath, refBed = args
    samplename = os.path.basename(bedFile).split('.')[0]
    print 'Running %s..' %samplename
    outFile = resultpath + '/' + samplename + '.tsv'
    bedString = open(refBed,'ru').readlines()
    refBed = BedTool(('\n'.join(bedString[1:])),from_string=True)
    coverageBed = BedTool(refBed).coverage(counts=True,b=bedFile)
    with open(outFile,'w') as ofile:
        header = bedString[0].strip() + '\t' + 'coverage\n'
        ofile.write(header)
        for line in coverageBed:
            ofile.write('\t'.join(line.fields).strip() + '\n')
        print 'Written %s' %ofile.name
    return 0

def main():
    projectpath = '/scratch/cdw2854/plasmaDNA'
    bedpath = projectpath + '/bedFiles'
    resultpath = projectpath + '/tissueContribution'
    os.system('mkdir -p %s' %resultpath)
    refpath = '/scratch/cdw2854/plasmaDNA/CTCFdata'
    refBed = refpath + '/cellCTCF.bed'
    bedFiles = glob.glob(bedpath + '/*.bed')
    bedFiles = filter(lambda x: re.search('NT|SRR|RNase|RNA-OCD|refA1',x),bedFiles)
    Pool(14).map(coverage,[(bedFile, resultpath, refBed) for bedFile in bedFiles])
    return 0


if __name__ == '__main__':
    main()
