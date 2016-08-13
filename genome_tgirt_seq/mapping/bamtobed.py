#!/stor/work/Lambowitz/cdw2854/src/miniconda2/bin/python

from pybedtools import BedTool
from pybedtools.bedtool import BAM
from pybedtools.cbedtools import Interval
import numpy as np
import fileinput
import sys

def filter_bed(bedline, min_length, max_length):
    chrom1,start1,end1 = bedline[:3]
    chrom2, start2, end2 = bedline[3:6]
    strand1, strand2 = bedline[8:10]
    name = bedline[6]
    start = np.min(map(int,[start1,start2]))
    end = np.max(map(int,[end1,end2]))
    length = end - start
    if strand2 != strand1 and min_length < length < max_length and chrom1==chrom2:
        return Interval(chrom = chrom1,
                start = start,
                end = end,
                name = name,
                score = length,
                strand = strand1)

def main():
    args = sys.argv
    if len(args)!=2:
        sys.exit('usage: python %s <bamfile> /dev/stdin for stdin' %args[0])
    min_length = 10
    max_length = 10000
    for frag in BedTool()\
            .bam_to_bed(i=args[1],bedpe=True,mate1=True, stream=True)\
            .each(filter_bed, min_length, max_length):
        print str(frag).strip()

if __name__=='__main__':
    main()
