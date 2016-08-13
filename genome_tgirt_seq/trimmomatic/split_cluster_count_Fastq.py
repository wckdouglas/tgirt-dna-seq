#!/usr/bin/env python

from Bio import SeqIO
import glob
import os
import gzip as gz
from itertools import izip
from multiprocessing import Pool

def split_files(arg):
    in_R1, cutoff, output_path = arg
    in_R2 = in_R1.replace('1P','2P')
    samplename = os.path.basename(in_R1).replace('1P.fq.gz','') + str(cutoff) + '_'
    out_r1 = output_path +'/' + samplename + '1P.fq.gz' 
    out_r2 = out_r1.replace('1P','2P')
    print 'Running %s' %samplename
    print 'Writing %s' %(out_r1)
    with gz.open(out_r1,'wb') as new_r1, gz.open(out_r2,'wb') as new_r2,\
	    gz.open(in_R1,'rb') as trimmed_r1, gz.open(in_R2,'rb') as trimmed_r2:
	counter = 0
	for r1, r2 in izip(SeqIO.parse(trimmed_r1,'fastq'), SeqIO.parse(trimmed_r2,'fastq')):
	    number_of_member = int(r1.description.split(' ')[1])
	    if number_of_member >= cutoff:
		new_r1.write(r1.format('fastq'))
		new_r2.write(r2.format('fastq'))
		counter += 1
    print 'For cluster member > %i, extracted %i sequences' %(cutoff,counter)

def main():
    project_path = '/scratch/02727/cdw2854/jurkatCells'
    trimmed_path= project_path +'/trimmed'
    output_path = project_path + '/filtered_clustered_fastq'
    if not os.path.isdir(output_path):
	os.mkdir(output_path)
    R1_files = glob.glob(trimmed_path + '/*1P.fq.gz')
    cutoff = range(0,11,2)
    args = [(R1, cut, output_path) for cut in cutoff for R1 in R1_files]
    pool = Pool(24)
    results = pool.map(split_files,args)
    pool.close()
    pool.join()
    return 0

if __name__ == '__main__':
    main()
