#!/usr/bin/env python

import subprocess
from multiprocessing import Pool
import glob
import os
import sys

projectpath = '/Users/wckdouglas/plasmaDNA'
resultpath = projectpath + '/results/plasmaResults'
bedpath = resultpath + '/bedFile' # fragment bed 6 for paried reads
result_bed = bedpath + '/noCTCF_bedFile'
refpath = projectpath + '/reference'
ctcf_bed = refpath + '/openChrom.gff'
genesBed = refpath + '/genes_count.bed'

def runCommand(command, samplename):
    shellCommand = subprocess.call(command, shell=True)
    message = '[%s] %s\n' %(samplename, command)
    sys.stderr.write(message)
    return 0


def countAndFilterBed(bedFile):
    samplename = os.path.basename(bedFile).split('.')[0]
    tempDir = result_bed + '/' + samplename
    runCommand('mkdir -p %s' %tempDir, samplename)
    command = "awk '$5 < 80' %s " %(bedFile) +\
	    "| intersectBed -a - -b %s -v " %ctcf_bed +\
            "| tee %s/%s_noCTCF.bed " %(result_bed, samplename) +\
            "| intersectBed -a - -b %s -wo " %(genesBed) +\
            "| awk '$6 == $12 && $NF > $5*0.5 {print $(NF-1), $(NF-2), $(NF-5) }'" +\
            "| sort -T %s " %tempDir +\
            '| uniq -c ' + \
            "| awk '{print $2,$3,$4,$1}' OFS='\\t'" +\
            '> %s/%s_count.tsv' %(result_bed, samplename)
    runCommand(command, samplename)
    runCommand('rm -rf %s' %tempDir, samplename)
    return 0

def main():
    bedFiles = glob.glob(bedpath + '/*bed')
    p = Pool(24).map(countAndFilterBed, bedFiles)
    return 0

if __name__ == '__main__':
    main()
