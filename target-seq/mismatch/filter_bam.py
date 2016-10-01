#!/usr/bin/env python

import pyximport
import pysam
pyximport.install(setup_args={
        'extra_link_args':pysam.get_libraries(),
        'include_dirs':pysam.get_include(),
        'define_macros':pysam.get_defines()
        })
from bam_filter import filterBAM
import sys

def pyfilterBAM(infile, outfile):
    in_bam = pysam.Samfile(infile,'rb')
    out_bam = pysam.Samfile(outfile,'w', template = in_bam)

    count = 0
    for aln in in_bam:
        if not aln.is_unmapped:
            if validateAlignment(aln):
                out_bam.write(aln)
                count += 1
    
    in_bam.close()
    out_bam.close()
    return count

arguments = sys.argv
if len(arguments) != 3:
    sys.exit('[usage]: python %s <in_bam> <out_bam>' %arguments[0])
infile = arguments[1]
outfile = arguments[2]
print 'Reading from %s' %infile
print 'Writing to %s' %outfile
count = filterBAM(infile, outfile)
print 'Written %i alignments' %(count)
