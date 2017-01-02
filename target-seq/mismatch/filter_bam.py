#!/usr/bin/env python

import pyximport
import pysam
import sys


def validate_alignment(aln):
    softclipped = 'S' in aln.cigarstring
    hardclipped = 'H' in aln.cigarstring
    indel = 'I' in aln.cigarstring or 'D' in aln.cigarstring
    mismatch_looks_good = aln.get_tag('NM') < 3
    on_hist1h3b = aln.reference_name == 'ENST00000621411'
    return  not softclipped and not hardclipped #and mismatch_looks_good #and not indel



def filter_bam(infile, outfile):
    in_bam = pysam.Samfile(infile,'rb')
    out_bam = pysam.Samfile(outfile,'wb', template = in_bam)

    count = 0
    for aln in in_bam:
        if not aln.is_unmapped:
            if validate_alignment(aln):
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
count = filter_bam(infile, outfile)
print 'Written %i alignments' %(count)
