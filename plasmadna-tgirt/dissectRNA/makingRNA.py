#!/bin/env python

from pybedtools import BedTool
import glob
import os
from multiprocessing import Pool

def len_filter(feature, lowerBound, upperBound):
    fragSize = feature.length
    return  fragSize > lowerBound and fragSize < upperBound

def mapFile(args):
    bedFile, motifBed, rnaBed, resultpath = args
    samplename = os.path.basename(bedFile).split('.')[0]
    print 'Running %s...' %samplename
    noDNABedname = '%s/%s.filtered.bed' %(resultpath, samplename)
    bed = BedTool(bedFile)
    bed.filter(len_filter, 30, 80) \
        .intersect(b = motifBed, v=True)\
        .intersect(b = rnaBed, f=0.6, r=True, wao=True) \
        .saveas(noDNABedname)
    print 'Written %s..' %noDNABedname
    
def makedir(directory):
    os.system('mkdir -p %s ' %directory)
    print 'Made directory: %s ' %directory
    return 0

def main():
    projectpath = '/scratch/cdw2854/plasmaDNA'
    bedpath = projectpath + '/dupbedFiles'
    resultpath = projectpath + '/smallReads'
    referencepath = '/scratch/cdw2854/reference'
    motifBed = referencepath + '/MotifFeatures.gff'
    rnaBed = referencepath + '/genes_count.bed'
    bindingSites = BedTool(motifBed)
    bedFiles = glob.glob(bedpath + '/*bed') 

    map(makedir, [resultpath])

    #Run bed files
    Pool(24).map(mapFile, [(bedFile, motifBed, rnaBed, resultpath) for bedFile in bedFiles])
    #map(mapFile, [(bedFile, motifBed, rnaBed, resultpath) for bedFile in bedFiles])


if __name__ == '__main__':
    main()
