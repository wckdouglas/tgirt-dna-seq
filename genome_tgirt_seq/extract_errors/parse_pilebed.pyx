#!/usr/bin/env python

import fileinput
from collections import defaultdict
from operator import itemgetter
import pandas as pd
import sys
import os

def summarize_mismatch(filename):
    cdef:
        str line, ref, A, C, T, G, read, mismatch_id
        int lineno, count

    read_ref_dict = defaultdict(int)
    bases = list('ACTG')
    samplename = os.path.basename(filename).split('.')[0]
    sys.stderr.write('Analyzing %s\n' %samplename)
    with open(filename, 'r') as infile:
        for lineno, line in enumerate(infile):
            if lineno != 0:
                fields = line.split('\t')
                ref, A, C, T, G = itemgetter(3,6,7,8,9)(fields)
                for read, count in zip(bases, map(int, [A,C,T,G])):
                    mismatch_id = ref+'>'+read
                    read_ref_dict[mismatch_id] += int(count)
                if lineno % 1000000 == 0:
                    sys.stderr.write('[%s] Parsed: %i lines\n' %(samplename, lineno))

    df = pd.DataFrame({
            'mismatch': read_ref_dict.keys(),
            'counts': read_ref_dict.values()}) \
        .assign(samplename = samplename)
    return df
