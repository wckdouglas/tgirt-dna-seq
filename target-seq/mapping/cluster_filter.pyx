import pysam
from pysam.calignmentfile cimport AlignmentFile, AlignedSegment
import os

cpdef int clusterCount(AlignedSegment aln):
    cdef:
        str name
        int cluster_count

    name = aln.query_name
    cluster_count = int(name.split('_')[-2])
    return cluster_count

cpdef int filterBam(str result_path, int min_read_cluster, str bam_file):
    cdef:
        str samplename
        str filtered_bam
        AlignmentFile in_bam, out_bam
        AlignedSegment aln

    samplename = bam_file.split('/')[-1].split('.')[0]
    filtered_bam = result_path + '/' + samplename + '.filtered.%i.bam' %(min_read_cluster)

    with pysam.Samfile(bam_file,'rb') as in_bam:
        with pysam.Samfile(filtered_bam,'wb',template = in_bam) as out_bam:
            [out_bam.write(aln) for aln in in_bam if clusterCount(aln) >= min_read_cluster]
    print 'Written %s' %filtered_bam
    os.system('samtools index %s' %filtered_bam)
    print 'Indexed %s' %filtered_bam
    return 0
