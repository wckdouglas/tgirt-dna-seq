#!/bin/env python

import argparse
import subprocess
import os
import sys

def getOpt():
    parser = argparse.ArgumentParser(description='Pipeline for trimming, mapping paired end plasma DNA')
    parser.add_argument('-1','--fq1',required=True, help='read1 fastqfile [string]')
    parser.add_argument('-o','--outdir',required=True, help='output directory')
    parser.add_argument('-x','--index',required=True, help='bwa index or hisat2 index')
    parser.add_argument('-a','--adaptor', default='adaptor.fa', help='Fasta file containing adaptor sequneces (default: adaptor.fa)')
    parser.add_argument('-t','--threads',default=1, type=int, help='Threads to be used (default=: 1)')
    args = parser.parse_args()
    return args

# running in shell
def runProcess(command, samplename):
    sys.stderr.write('[%s] %s\n' %(samplename, command))
    result = subprocess.call('time ' + command, shell=True)
    return 0

#Trimming
def trimming(fq1, threads, trim_path, samplename, adaptor):
    sys.stderr.write('Running trim process with %s\n' %samplename)
    #ILLUMINACLIP:<fastaWithAdaptersEtc>:<seed mismatches>:<palindrome clip
    #       threshold>:<simple clip threshold>:<minAdapterLength>:<keepBothReads>
    options='ILLUMINACLIP:%s:2:10:10:1:true ' %(adaptor)+\
	    'LEADING:10 TRAILING:10 '  +\
	    'SLIDINGWINDOW:4:8  MINLEN:18  AVGQUAL:15'
    command = 'time trimmomatic ' +\
        'PE -threads %i '  %(threads)+\
        '-basein %s ' %(fq1) + \
        '-baseout %s/%s.fq.gz ' %(trim_path, samplename) + \
        options
    #runProcess(command,samplename)
    return 0

#MAPPING
def mappingProcess(samplename, trim_path, index, threads, bam_path):
    sys.stderr.write('Running mapping with %s\n' %samplename)
    file1 = trim_path + '/' + samplename + '_1P.fq.gz'
    file2 = file1.replace('1P','2P')
    bam_file = '%s/%s.bam' %(bam_path, samplename)
    command = 'bowtie2 --threads %i --no-discordant  --local ' %(threads)+\
	    '--no-mixed -x %s -1 %s -2 %s ' %(index, file1, file2 ) +\
            '| samtools view -@ %i -b ' %(threads) +\
            '> %s' %bam_file
    runProcess(command, samplename)
    return bam_file

def makeBed(bam_file, samplename, bed_path):
    sys.stderr.write('Running post mapping processes with %s\n' %samplename)
    bed_file = '%s/%s.bed' %(bed_path, samplename)
    command = 'bamtools filter -script flag_filter.json -in %s' %(bam_file)+\
	'| samtools fixmate -r - -' + \
        '| bedtools bamtobed -mate1 -bedpe '+\
        '| python bedpetobed.py - '+\
        '> %s' %(bed_file)
    runProcess(command,samplename)
    return bed_file

def makeRmdupBed(bed_file, samplename, rmdup_bed_path):
    sys.stderr.write('Running post mapping processes with %s\n' %samplename)
    temp_dir = rmdup_bed_path+ '/' + samplename
    makedir(temp_dir)
    command = 'cat %s ' %(bed_file) +\
        '| sort -k1,1 -k2,2n -k3,3n -k6,6  --temporary-directory=%s -u ' %(temp_dir)+\
        '> %s/%s.bed' %(rmdup_bed_path, samplename)
    runProcess(command,samplename)
    os.rmdir(temp_dir)
    return 0

def makedir(directory):
    if not os.path.isdir(directory):
        os.mkdir(directory)
        sys.stderr.write('Make directory %s\n' %directory)
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
    bed_path = outdir + '/bedFiles'
    rmdup_bed_path = outdir + '/rmdup_bed'
    map(makedir,[trim_path, bam_path, rmdup_bed_path, bed_path])

    #trim
    trim = trimming(fq1, threads, trim_path, samplename, adaptor)
    #map
    bam_file = mappingProcess(samplename, trim_path, index, threads, bam_path)
    bed_file = makeBed(bam_file, samplename, bed_path)
    makeRmdupBed(bed_file, samplename, rmdup_bed_path)
    sys.stderr.write('Finished mapping %s\n' %samplename)
    return 0

if __name__ == '__main__':
    args = getOpt()
    main(args)
