#!/usr/bin/env python

import pandas as pd
import numpy as np


ref_path = '/stor/work/Lambowitz/ref/hg19/methylation'
ref_table = 'pnas.1508736112.sd01.xlsx'
tissue_methylation = pd.read_excel(ref_path + '/' + ref_table, 
                                   sheet = 0) 
coor_names = ['chrom','start','end','Genomic location','Tissue']
tissue_methylation = tissue_methylation \
        .assign(chrom = lambda d: map(lambda x: x.split(':')[0], d['Genomic location']))\
        .assign(start = lambda d: map(lambda x: x.split(':')[1].split('-')[0], d['Genomic location'])) \
        .assign(end = lambda d: map(lambda x: x.split('-')[1], d['Genomic location'])) \
        .pipe(lambda d: d[coor_names])
out_bed_name = ref_path + '/methyl_table.bed'
tissue_methylation.to_csv(out_bed_name,sep='\t', 
                          header=False, index=False)
print 'Written: ', out_bed_name

