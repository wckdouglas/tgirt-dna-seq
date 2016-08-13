#!/usr/bin/env python

import numpy as np
from scipy import fftpack
import seaborn as sns

def findIntercepts(wpsArray, outFile, chromosome,samplename):
    '''
    looking for the positions that across 0 in the de-noised signal
    '''
    wpsArray = np.asarray(wpsArray)

    signs = np.sign(wpsArray)
    signs[signs==0] = -1
    start = np.where(np.diff(signs)>0)[0]  
    end = np.where(np.diff(signs)<0)[0]  
    numberOfPeak = 0
    lowPeakBound, upperPeakBound = (15, 120)
    for i, j  in zip(start, end):
        nucleosomeSize = np.abs(j - i)
        coverageScore = np.max(wpsArray[i:j])
        if nucleosomeSize > lowPeakBound and nucleosomeSize < upperPeakBound:
            numberOfPeak += 1
            peakname = '%s_peak_%i' %(chromosome,numberOfPeak)
            peakPos = (i + j)/2
            line = '\t'.join(map(str,[chromosome, i, j, peakname, peakPos ,'+', coverageScore]))
            outFile.write(line+'\n')
    print 'Written %i peaks to %s' %(numberOfPeak, outFile.name)
    return 0

def filterWPS(wps):
    intensity = fftpack.fft(wps)
    freq = fftpack.fftfreq(len(wps))
    periodicity = 1/(freq)
    intensity[(periodicity > 120) & (periodicity < 200)] = 0
    intensity[(periodicity < -120) & (periodicity > -200)] = 0
    filtered = fftpack.ifft(intensity)
    return filtered

def readWig(args):
    wigFile, resultpath = args 
    basename = os.path.basename(wigFile)
    samplename = basename.split('.')[0]
    chrom = basename.split('.')[1]
    print 'Running: %s ' %basename
    wps = np.array([float(line.split('\t')[1]) for line in open(wigFile,'ru') if 'chrom' not in line],
            dtype=np.float64)
    filteredWPS = filteredWPS(wps)
    with open('%s/%s.%s.bed' %(resultpath, samplename, chrom)) as outFile:
        print 'Writing %s' %(outFile.name)
        findIntercepts(filteredWPS, outFile, chrom, samplename)
    return 0

def main():
    projectpath = '/scratch/02727/cdw2854/plasmaDNA'
    datapath = projectpath + '/genomeWPS'
    resultpath = projectpath + '/rnaBed'
    if not os.path.isdir(resultpath):
        os.mkdir(resultpath)
    wigFiles = glob.glob(datapath + '/*merged.wig') # wig file storing wps scores for each genome position
    map(readWig, [(wigFile, resultpath) for wigFile in wigFiles])

if __name__ == '__main__':
    main()
