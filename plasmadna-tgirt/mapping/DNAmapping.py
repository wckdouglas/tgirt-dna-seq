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
    options='ILLUMINACLIP:%s:2:10:10:1:true MINLEN:15' %(adaptor)
    command = 'java -Xmx1G -jar trimmomatic-0.36.jar ' +\
        'PE -threads %i '  %(threads)+\
        '-basein %s ' %(fq1) + \
        '-baseout %s/%s.fq.gz ' %(trim_path, samplename) + \
        options
    runProcess(command,samplename)
    return 0

#MAPPING
def bwa_mapping(samplename, trim_path, index, threads, bam_path):
    '''
    map reads with bwa mem
    '''
    sys.stderr.write('Running mapping with %s\n' %samplename)
    file1 = trim_path + '/' + samplename + '_1P.fq.gz'
    file2 = file1.replace('1P','2P')
    bam_file = '%s/%s.bam' %(bam_path, samplename)
    #runProcess('bwa shm %s' %(index), samplename)
    command = 'bwa mem -t %i ' %(threads)+\
	    '%s %s %s ' %(index, file1, file2 ) +\
            '| samtools view -@ %i -b ' %(threads) +\
            '> %s' %bam_file
    runProcess(command, samplename)
    return bam_file

def make_bed(bam_file, samplename, bed_path):
    '''
    merge pair end reads into fragments and convert to bed file
    '''
    sys.stderr.write('Running post mapping processes with %s\n' %samplename)
    bed_file = '%s/%s.bed' %(bed_path, samplename)
    #command = 'bamtools filter -script flag_filter.json -in %s' %(bam_file)+\
    command = 'samtools view -bF2048 -F1024 -F512 -F256 -F4 -F8 %s' %(bam_file)+\
	'| samtools fixmate -r - -' + \
        '| bedtools bamtobed -mate1 -bedpe '+\
        "| awk '$1!=\".\" && $NF!=$(NF-1) && $1==$4 "+\
                "{start=$2;end=$3} {if($5<start) start=$5} {if($6>end) end=$6} "+\
                "{print $1,start,end,$7,end-start,$9}' OFS='\\t'"+\
        "| awk '$(NF-1) < 1000'" +\
        '> %s' %(bed_file)
    runProcess(command,samplename)
    return bed_file

def split_bed(bed_file, samplename, split_bed_path):
    '''
    Splitting bed by chromosome
    '''
    chroms = range(1,22)
    chroms.extend(list('XY'))
    chroms.append('MT')
    for chrom in chroms:
        out_name = '%s/%s.%s.bed' %(split_bed_path, samplename, chrom)
        if os.path.isfile(out_name):
            os.remove(out_name)
    command = 'awk \'$1~/^[0-9]+$|^[XY]$|^MT$/ {print $0 >> "%s/%s."$1".bed"}\' %s' \
            %(split_bed_path, samplename,bed_file)
    runProcess(command, samplename)
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
    samplename = os.path.basename(fq1).replace('_R1','').replace('_001','').split('.')[0]

    #makedir
    trim_path= outdir + '/trimmedFiles'
    bam_path= outdir + '/bamFiles'
    bed_path = outdir + '/bedFiles'
    rmdup_bed_path = outdir + '/rmdup_bed'
    split_bed_path = bed_path + '/chrom_split_bed'
    map(makedir,[trim_path, bam_path, rmdup_bed_path, bed_path, split_bed_path])

    #trim
    trim = trimming(fq1, threads, trim_path, samplename, adaptor)
    #map
    bam_file = bwa_mapping(samplename, trim_path, index, threads, bam_path)
    bed_file = make_bed(bam_file, samplename, bed_path)
    split_bed(bed_file, samplename, split_bed_path)
    sys.stderr.write('Finished mapping %s\n' %samplename)
    return 0

if __name__ == '__main__':
    args = getOpt()
    main(args)
