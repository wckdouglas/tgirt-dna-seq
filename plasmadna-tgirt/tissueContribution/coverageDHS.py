#!/usr/bin/env python

from pybedtools import BedTool, set_tempdir
import os
import glob
from multiprocessing import Pool
import re
from functools import partial

def coverage(resultpath, ref_bed_name, header, bedFile):
    samplename = os.path.basename(bedFile).split('.')[0]
    print 'Running %s..' %samplename
    outFile = resultpath + '/' + samplename + '.tsv'
    sorted_bed =  BedTool(bedFile).sort()
    with open(outFile,'w') as out_file:
        out_file.write(header)
        for line in BedTool(ref_bed_name)\
                .intersect(b = sorted_bed.fn, wb=True, f=0.3, r=True, sorted=True):
            ref_line = line.fields[:24]
            score = line.fields[-3]
            ref_line.append(score)
            out_file.write('\t'.join(ref_line) + '\n')
        print 'Written %s' %out_file.name
    return 0

def main():
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bedpath = projectpath + '/genomeWPS/bed_files'
    resultpath = projectpath + '/tissueContribution'
    os.system('mkdir -p %s' %resultpath)
    set_tempdir(resultpath)
    refpath = os.environ['REF'] + '/ctcfData'
    refBed = refpath + '/cellCTCF.bed'
    bedFiles = glob.glob(bedpath + '/*.Short.bed')
    bedFiles = filter(lambda x: re.search('SRR|P1022',x),bedFiles)

    #process reference bed
    bedString = open(refBed,'ru').readlines()
    ref_bed_name = resultpath + '/dhs_ref.bed'
    BedTool(('\n'.join(bedString[1:])),from_string=True)\
        .sort()\
        .moveto(ref_bed_name)
    header = bedString[0].strip() + '\t' + 'WPS\n'

    coverage_func = partial(coverage, resultpath, ref_bed_name, header)
    Pool(12).map(coverage_func, bedFiles)
    #map(coverage_func, bedFiles)
    return 0


if __name__ == '__main__':
    main()
