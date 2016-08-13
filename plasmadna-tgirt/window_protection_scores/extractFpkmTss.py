#!/bin/env python

import fileinput

class transcript:
    def __init__(self,chrom,start,FPKM,strand,id,name):
        self.chrom = chrom
        self.start = int(start)
        self.exprs = FPKM
        self.strand = strand
        self.id = id
        self.name = name

for line in fileinput.input():
    if line[0] !='#':
        line = line.strip()
        columns = line.split('\t')
        if columns[2] == 'transcript':
            info = columns[8]
            info = ''.join(info.split(' '))
            transcriptInfo = {}
            for i in info.split(';') :
                idx = i.split('"')
                if len(idx)>1:
                    transcriptInfo[idx[0]] = idx[1]
        trans = transcript(columns[0],columns[3],transcriptInfo['FPKM'],columns[6],transcriptInfo['reference_id'],transcriptInfo['ref_gene_id'])
        print '\t'.join([trans.chrom,str(trans.start-1000),str(trans.start+1000),trans.id,'0',trans.strand,trans.name,trans.exprs])
