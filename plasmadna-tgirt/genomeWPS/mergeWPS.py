#!/usr/bin/env python -u

import glob
import os
from multiprocessing import Pool
from itertools import izip

def mergeFile(args):
    shortFile, longFile, resultpath = args
    samplename = os.path.basename(longFile).split('.')[0]
    chrom = os.path.basename(longFile).split('.')[1]
    outfile = resultpath + '/' + samplename + '.' + chrom + '.merged.wig'
    print 'Running sample: %s for chrom %s ' %(samplename, chrom)
    lineno = 0
    with open(outfile,'w') as of, open(shortFile,'r') as sf, open(longFile,'r') as lf:
        for shortLine, longLine in izip(sf, lf):
            if lineno == 0: 
                assert shortLine == longLine, 'Wrong chromosome!! %s' %samplename
                of.write(shortLine)
            else:
                shortFields = shortLine.strip().split('\t')
                longFields = longLine.strip().split('\t')
                shortWPS = float(shortFields[1])
                longWPS = float(shortFields[1])
                shortPos = shortFields[0]
                longPos = longFields[0]
                assert longPos == shortPos, 'Wrong position!! %s' %samplename
                line = '%s\t%.3f\n' %(longPos, shortWPS + longWPS)
                of.write(line)
            lineno += 1
    print 'Finished: %s for chrom %s' %(samplename, chrom)
    return 0

def main():
    projectpath = '/scratch/02727/cdw2854/plasmaDNA/'
    datapath = projectpath + 'genomeWPS'
    resultpath = projectpath + 'genomeWPS'
    longWpsFiles = glob.glob(datapath + '/*Long.wig')
    shortWpsFiles = [longWpsFile.replace('.Long.','.Short.') for longWpsFile in longWpsFiles]
    pool = Pool(processes=48, maxtasksperchild=1)
    pool.map_async(mergeFile, [(swf, lwf, resultpath) for swf, lwf in zip(shortWpsFiles,longWpsFiles)]).get()
    return 0

if __name__ == '__main__':
    main()
