#!/usr/bin/env python

import fileinput
import sys

read = ''

def clip_seq(line, read):
    if read == 'read1':
        print line[1:]
    elif read == 'read2':
        print line[4:]



for i, line in enumerate(fileinput.input()):
    line = line.strip()
    if line[0] == '@':
        line_count = 0
        if line.endswith('/1'):
            read = 'read1'
        elif line.endswith('/2'):
            read = 'read2'
        print line
    else:
        line_count += 1
        if line_count == 1:
            clip_seq(line, read)
        elif line_count == 2:
            if line == '+':
                print '+'
            else:
                sys.exit('Wrong Parse')
        elif line_count == 3:
            clip_seq(line, read)
