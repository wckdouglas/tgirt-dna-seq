#!/bin/env python

import argparse
import subprocess
import os
import sys

def getOpt():
    parser = argparse.ArgumentParser(description='Pipeline for trimming, mapping paired end plasma DNA')
    parser.add_argument('-1','--fq1',required=True, help='read1 fastqfile [string]')
    parser.add_argument('-o','--outdir',required=True, help='output directory')
    parser.add_argument('-x','--index',required=True, help='bwa ')
    parser.add_argument('-a','--adaptor', default='adaptor.fa', help='Fasta file containing adaptor sequneces (default: adaptor.fa)')
    parser.add_argument('-t','--threads',default=1, type=int, help='Threads to be used (default=: 1)')
    args = parser.parse_args()
    return args

# running in shell
def runProcess(command, samplename):
    sys.stderr.write('[%s] Running %s\n' %(samplename, command))
    result = subprocess.call('time ' + command, shell=True)
    return 0

#Trimming 
def trimming(fq1, threads, trim_path, samplename, adaptor):
    sys.stderr.write('Running trim process with %s\n' %samplename)
    options='ILLUMINACLIP:%s:2:10:10:1:true ' %(adaptor)+\
	    'LEADING:10 TRAILING:10 SLIDINGWINDOW:4:8  MINLEN:18  AVGQUAL:26'
    command = 'trimmomatic PE -threads %i '  %(threads)+\
        '-basein %s ' %(fq1) + \
        '-baseout %s/%s.fq.gz ' %(trim_path, samplename) + \
        options
    runProcess(command,samplename)
    return 0

#MAPPING
def mappingProcess(samplename, trim_path, index, threads, rmdup_bam_path, bam_path):
    file1 = trim_path + '/' + samplename + '_1P.fq.gz'
    file2 = file1.replace('1P','2P')
    postMapProcess = '| samtools view -b@ %i ' %(threads) +\
                '> %s/%s.bam' %(bam_path, samplename)
    command = 'bwameth.py  --threads=%i ' %(threads)+\
		'--reference=%s %s %s ' %(index, file1, file2 )  +\
                '| samtools view -@ %i -b ' %(threads) +\
                '| samtools sort -@ %i -T %s/%s -O bam' %(threads, bam_path, samplename) +\
                '> %s/%s.bam ' %(bam_path,samplename) 
    runProcess(command, samplename)
    return '%s/%s.bam' %(bam_path, samplename)

def postProcess(bamFile, samplename, rmdup_bam_path):
    command = 'picard MarkDuplicates ' + \
            'ASSUME_SORT_ORDER=coordinate ' +\
	    'INPUT=%s ' %(bamFile) +\
    	    'METRICS_FILE=%s/%s.txt ' %(rmdup_bam_path, samplename) +\
	    'OUTPUT=%s/%s.bam ' %(rmdup_bam_path, samplename)
    runProcess(command,samplename)
    return 0

def makedir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
        sys.stderr.write('Make direcroty %s\n' %directory)
    return 0

def main(args):
    fq1 = args.fq1
    outdir = args.outdir
    index = args.index
    adaptor = args.adaptor
    threads = args.threads

    # set up variables
    suffix = '.'.join(fq1.split('.')[-2:]) if fq1.split('.')[-1] == 'gz' else fq1.split('.')[-1]
    samplename = os.path.basename(fq1).replace('_R1_001','').split('.')[0]

    #makedir
    trim_path= outdir + '/trimmedFiles'
    bam_path= outdir + '/bamFiles'
    bed_path = outdir + 'bedFiles'
    rmdup_bam_path = outdir + '/rmdupBAM'
    map(makedir,[trim_path, bam_path, rmdup_bam_path])

    ####triming sequences
    #trim = trimming(fq1, threads, trim_path, samplename,  adaptor)
    #map
    bamFile = mappingProcess(samplename, trim_path, index, threads, rmdup_bam_path, bam_path)
    ok = postProcess(bamFile, samplename, rmdup_bam_path)
    return 0

if __name__ == '__main__':
    args = getOpt()
    main(args)
