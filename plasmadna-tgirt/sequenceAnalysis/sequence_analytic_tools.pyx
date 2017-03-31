from collections import defaultdict
from operator import itemgetter
from pyfaidx import Fasta

def make_dinucleotide(sequence):
    '''
    output dinucleotides from a sequence
    '''
    cdef:
        str first_nucleotide = sequence[:-1]
        str second_nucleotide = sequence[1:]
        str x,y

    for x, y in zip(first_nucleotide, second_nucleotide):
        yield x+y


def extract_dinucleotide(dinucleotide_count_dict, sequence):
    '''
    for each posiiton add dinucleotide composition
    '''
    cdef:
        int i
        str di

    dinucleotides = make_dinucleotide(sequence)
    for i, di in enumerate(dinucleotides):
        dinucleotide_count_dict[di][i] += 1

    return dinucleotide_count_dict

def fragment_center(start, end, window):
    '''
    Given the fragment ends, calculate center and output window
    '''
    cdef:
        long center, center_start, center_end

    center = (long(end) + long(start))
    center_start = center - window
    center_end = center + window
    return center_start, center_end

def parse_bed(bed_file, ref_fasta, window, regular_chrom):
    '''
    extract dinucleotide for center of 167-nt fragments
    '''
    cdef:
        int fragment_counts = 0
        str chrome, strand, fragment

    dinucleotide_count_dict = defaultdict(lambda : defaultdict(int))
    ref = Fasta(ref_fasta)
    with open(bed_file, 'r') as in_bed:
        for fragment in in_bed:
            fields = fragment.rstrip().split('\t')
            chrom, start, end, insert_size, strand = itemgetter(0,1,2,4,5)(fields)
            if int(insert_size) == 167 and chrom in regular_chrom:
                center_start, center_end = fragment_center(start, end, window)
                if center_end < len(ref[chrom]):
                    fragment_counts += 1
                    sequence = ref[chrom][center_start:center_end]
                    sequence = sequence if strand == '+' else sequence.reverse.complement
                    dinucleotides_count_dict = extract_dinucleotide(dinucleotide_count_dict, str(sequence))
    print 'Parsed %i 167-nt fragments' %(fragment_counts)
    return dinucleotide_count_dict
