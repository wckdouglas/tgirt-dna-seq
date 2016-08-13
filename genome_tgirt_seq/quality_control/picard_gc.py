#!/bin/env python

import glob
import os
from functools import partial
import time

def runProcess(command):
    print command
    os.system(command)

def outBamName(result_path, samplename, analysis_type):
    return result_path + '/' + samplename + '.' + analysis_type + '.bam'

def filterBam(in_bam, outNameFunc):
    out_bam = outNameFunc('filter')
    command = 'bamtools filter -script flag_filter.json -in %s > %s' %(in_bam,out_bam)
    runProcess(command)
    return out_bam

def fixMate(in_bam, outNameFunc):
    out_bam = outNameFunc('fixedmate')
    command = 'picard FixMateInformation ' + \
        'INPUT=%s ' %(in_bam) +\
        'OUTPUT=%s '%(out_bam) +\
        'SORT_ORDER=coordinate ' +\
	'ADD_MATE_CIGAR=true '
    runProcess(command)
    return out_bam

def sortBam(in_bam, result_path, samplename):
    out_bam = '%s/%s.sorted.bam'  %(result_path, samplename)
    command = 'samtools sort -@ 12 -T %s/%s -O bam %s ' %(result_path, samplename, in_bam) +\
            '>  %s' %(out_bam)
    runProcess(command)
    return out_bam

def gcCollect(in_bam, figures_path, samplename, result_path, ref):
    command =' picard CollectGcBiasMetrics '+\
	'SCAN_WINDOW_SIZE=100 ' + \
	'INPUT=%s ' %(in_bam) +\
        'CHART_OUTPUT=%s/%s.pdf '  %(figures_path, samplename) + \
        'OUTPUT=%s/%s.txt ' %(result_path, samplename) +\
	'SUMMARY_OUTPUT=%s/%s.summary ' %(result_path, samplename) +\
	'REFERENCE_SEQUENCE=%s ' %ref +\
	'ASSUME_SORTED=true '
    runProcess(command)


def pipeline(result_path, figures_path, ref, bam_file):
    start = time.time()
    samplename = os.path.basename(bam_file).replace('.bam','')
#    outNameFunc = partial(outBamName,result_path, samplename)
#    filtered_bam = filterBam(bam_file, outNameFunc)
#    fixmate_bam = fixMate(filtered_bam, outNameFunc)
    sorted_bam = sortBam(bam_file, result_path, samplename)
    gcCollect(sorted_bam, figures_path, samplename, result_path, ref)
    end = time.time()
    time_lapsed = (end - start)/float(60)
    print 'Finished %s in %.3f min' %(samplename, time_lapsed)
    return 0

def main():
    project_path = '/stor/scratch/Lambowitz/cdw2854/genomeDNA'
    bam_path= project_path + '/bamFiles'
    result_path = project_path + '/picard_results'
    figures_path = project_path + '/figures'
    ref_path = '/stor/work/Lambowitz/ref/hg19/Sequence/WholeGenomeFasta'
    ref = ref_path + '/genome.fa'
    for path in [result_path, figures_path]:
        if not os.path.isdir(path):
            os.makedirs(path)

    bam_files = glob.glob(bam_path + '/*bam')
    picard_func = partial(pipeline, result_path, figures_path, ref)
    map(picard_func, bam_files)

if __name__ == '__main__':
    main()
