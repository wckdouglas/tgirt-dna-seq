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

arguments = sys.argv
if len(arguments) != 3:
    sys.exit('[usage]: python %s <in_bam> <out_bam>' %arguments[0])
infile = arguments[1]
outfile = arguments[2]
count = filterBAM(infile, outfile)
print 'Written %i alignments' %(count)
