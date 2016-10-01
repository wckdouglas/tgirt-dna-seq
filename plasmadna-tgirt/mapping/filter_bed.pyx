from pybedtools.cbedtools cimport Interval
from pybedtools.cbedtools import Interval
import sys
from cpython cimport bool

cpdef Interval filterBed(Interval bedline, int min_length, int max_length):

    cdef:
        unicode chrom1, chrom2
        unicode start1, start2
        unicode end1, end2
        unicode strand1, strand2
        unicode name
        long start, end, length
        bool strandeness_correct, length_correct, chrom_correct

    chrom1,start1,end1 = bedline[:3]
    chrom2, start2, end2 = bedline[3:6]
    strand1, strand2 = bedline[8:10]
    name = bedline[6]
    start = min(map(int,[start1,start2]))
    end = max(map(int,[end1,end2]))
    length = end - start

    strandeness_correct = strand1 != strand2
    length_correct = min_length < length < max_length
    chrom_correct = chrom1 == chrom2

    if strandeness_correct and length_correct and chrom_correct:
        alignment = Interval(chrom = chrom1,
                start = start,
                end = end,
                name = name,
                score = length,
                strand = strand1)
        return alignment


cpdef int processFile(bed_iterator):
    cdef:
        int min_length = 10
        int max_length = 10000
        Interval frag

    for frag in bed_iterator \
        .each(filterBed, min_length, max_length):
        print str(frag).strip()
    return 0
