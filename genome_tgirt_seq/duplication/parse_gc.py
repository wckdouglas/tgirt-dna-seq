#!/usr/bin/env python

import fileinput

for line in fileinput.input():
    count = '1'
    member_count, length, gc_percentage = line.rstrip().split('\t')
    if '_member' in member_count:
        count = member_count.split('_')[1]
    print '%s\t%s\t%s' %(count, length, gc_percentage)
