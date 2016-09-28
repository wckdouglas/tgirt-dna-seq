#!/usr/bin/env python

import sys
from pybedtools import BedTool
import pyximportcpp
pyximportcpp.install(setup_args={
            'include_dirs': ['include/']
            })
from filter_bed import filterBed, processFile

def main():
    args = sys.argv
    if len(args)!=2:
        sys.exit('usage: python %s <bedpe_file> /dev/stdin for stdin' %args[0])
    bed_file = args[1]
    handle = bed_file if bed_file != '/dev/stdin' and bed_file != '-' else sys.stdin
    bed_iterator = BedTool(handle)
    done = processFile(bed_iterator)


if __name__=='__main__':
    main()
