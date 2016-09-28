import pysam
from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
from cpython cimport bool


cpdef bool validateAlignment(AlignedSegment aln):
    cdef:
        bool softclipped, mismatch_looks_good, on_hist1h3b

    softclipped = 'S' in aln.cigarstring
    mismatch_looks_good = aln.get_tag('NM') < 4
    on_hist1h3b = aln.reference_name == 'ENST00000621411'
    return  not softclipped and mismatch_looks_good and on_hist1h3b


cpdef int filterBAM(infile, outfile):
    cdef:
        AlignmentFile in_bam
        AlignmentFile out_bam
        AlignedSegment aln
        int count

    in_bam = pysam.Samfile(infile,'rb')
    out_bam = pysam.Samfile(outfile,'wb', template = in_bam)

    count = 0
    for aln in in_bam:
        if not aln.is_unmapped:
            if validateAlignment(aln):
                out_bam.write(aln)
                count += 1

    
    in_bam.close()
    out_bam.close()
    return count
