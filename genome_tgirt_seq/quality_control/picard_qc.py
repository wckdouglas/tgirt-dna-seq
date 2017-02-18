#!/bin/env python

import glob
import os
from functools import partial
import time
from multiprocessing import Pool

def runProcess(command):
    print command
    os.system(command)

def outBamName(result_path, samplename, analysis_type):
    return result_path + '/' + samplename + '.' + analysis_type + '.bam'

def filterBam(in_bam, outNameFunc):
    out_bam = outNameFunc('filter')
    command = 'bamtools filter -script flag_filter.json ' +\
                '-in %s ' %(in_bam)+\
            '| samtools sort -O bam -T %s ' %out_bam.replace('.bam','')+\
            '> %s' %(out_bam)
    command2 = 'samtools index %s' %(out_bam)
    runProcess(command)
    runProcess(command2)
    return out_bam

def MarkDuplicates(in_bam, outNameFunc):
    out_bam = outNameFunc('MarkDuplicate')
    command = 'picard MarkDuplicates ' + \
        'INPUT=%s ' %(in_bam) +\
        'OUTPUT=%s '%(out_bam) +\
        'ASSUME_SORT_ORDER=coordinate ' +\
        'REMOVE_SEQUENCING_DUPLICATES=false ' +\
        'CREATE_INDEX=true '+\
	'METRICS_FILE=%s ' %(out_bam.replace('bam','duplicate.metric'))
    runProcess(command)
    return out_bam

def collect_gc(in_bam, figures_path, samplename, result_path, ref):
    command = 'picard CollectGcBiasMetrics '+\
	'SCAN_WINDOW_SIZE=100 ' + \
	'INPUT=%s ' %(in_bam) +\
        'CHART_OUTPUT=%s/%s.pdf '  %(figures_path, samplename) + \
        'OUTPUT=%s/%s.gc_metrics ' %(result_path, samplename) +\
	'SUMMARY_OUTPUT=%s/%s.summary ' %(result_path, samplename) +\
	'REFERENCE_SEQUENCE=%s ' %ref +\
	'ASSUME_SORTED=true ' +\
        'VALIDATION_STRINGENCY=SILENT'
    runProcess(command)

def subsampling(bam_file):
    fold = 500000
    subsampled_bams = []
    for seed in xrange(10):
        subsampled_bam = bam_file.replace('.bam','.%i.subsampled.bam' %(seed))
        command = 'samtools view -bF 256 -F 4 -F 1024 -F 2048 %s ' %(bam_file) +\
            '| bedtools sample -n %i -seed %i -i - ' %(fold, seed) +\
            '> %s ' %(subsampled_bam)
        runProcess(command)
        subsampled_bams.append(subsampled_bam)
    return subsampled_bams

def collect_wgs(bam_files, ref):
    for bam_file in bam_files:
        out_metric = bam_file.replace('.bam','.wgs_metrics')
        command = 'picard CollectWgsMetrics '+\
            'INPUT=%s ' %bam_file +\
            'COVERAGE_CAP=300 ' +\
            'COUNT_UNPAIRED=true ' +\
            'MINIMUM_BASE_QUALITY=0 '+\
            'MINIMUM_MAPPING_QUALITY=0 '+\
            'REFERENCE_SEQUENCE=%s ' %(ref)+\
            'OUTPUT=%s ' %(out_metric) +\
            'INCLUDE_BQ_HISTOGRAM=true'
        runProcess(command)

def collect_alignment(bam_files, ref):
    for bam_file in bam_files:
        out_metric = bam_file.replace('.bam','.alignment_metrics')
        command = 'picard CollectAlignmentSummaryMetrics ' +\
            'REFERENCE_SEQUENCE=%s ' %ref +\
            'INPUT=%s ' %bam_file +\
            'OUTPUT=%s ' %out_metric +\
            'ASSUME_SORTED=true ' 
        runProcess(command)

def pipeline(result_path, figures_path, ref, bam_file):
    start = time.time()
    samplename = os.path.basename(bam_file).replace('.bam','')
    outNameFunc = partial(outBamName,result_path, samplename)
    filtered_bam = filterBam(bam_file, outNameFunc)
    dedup_bam = MarkDuplicates(filtered_bam, outNameFunc)
    subsampled_bams = subsampling(dedup_bam)
    collect_gc(dedup_bam, figures_path, samplename, result_path, ref)
    collect_wgs(subsampled_bams, ref)
    collect_alignment(subsampled_bams, ref)
    

    end = time.time()
    time_lapsed = (end - start)/float(60)
    print 'Finished %s in %.3f min' %(samplename, time_lapsed)
    return 0

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome'
    bam_path= project_path + '/bamFiles'
    #project_path = '/stor/work/Lambowitz/cdw2854/genomeDNA'
    #bam_path = '/stor/work/Lambowitz/cdw2854/genomeDNA/bamFiles'
    result_path = project_path + '/picard_results'
    figures_path = project_path + '/figures'
    ref_path = '/stor/work/Lambowitz/ref/Ecoli'
    ref = ref_path + '/k12_mg1655.fa'
    #ref = '/stor/work/Lambowitz/ref/hg19/Sequence/all_seq/reference.fasta'
    #ref = '/stor/work/Lambowitz/ref/hg19/Sequence/WholeGenomeFasta/genome.fa'
    for path in [result_path, figures_path]:
        if not os.path.isdir(path):
            os.makedirs(path)

    bam_files = glob.glob(bam_path + '/*.bam')
    bam_files = filter(lambda x: 'pb' not in x, bam_files)
    picard_func = partial(pipeline, result_path, figures_path, ref)
    p = Pool(20)
    p.map(picard_func, bam_files)
    p.close()
    p.join()

if __name__ == '__main__':
    main()
