#!/usr/bin/env python

import pandas as pd
import numpy as np

def process_marker(ref_table, marker_type):
    marker_table = pd.read_excel(ref_table, 
                        sheetname='Type.%s.biomarkers' %(marker_type)) 
    tissue_methylation = marker_table \
        .assign(chrom = lambda d: map(lambda x: x.split(':')[0], d['Genomic location']))\
        .assign(start = lambda d: map(lambda x: x.split(':')[1].split('-')[0], d['Genomic location'])) \
        .assign(end = lambda d: map(lambda x: x.split('-')[1], d['Genomic location'])) 

    if marker_type == 'II':
        tissue_methylation['Tissue'] = 'none'
    return tissue_methylation

def main():
    ref_path = '/stor/work/Lambowitz/ref/hg19/methylation'
    ref_table = ref_path + '/pnas.1508736112.sd01.xlsx'
    marker_types = ['I','II']
    df = [process_marker(ref_table, marker_type )for marker_type in marker_types]
    df = pd.concat(df, axis=0)

    #write deconcolute table
    df.drop(['C.V.','S.D.','Placenta-informative','Z score','Maximum','Minimum'], 
            axis=1, inplace=True)
    tablename = ref_path + '/methyl_tissue_table.tsv'
    df\
        .drop(['chrom','start','end'], axis=1)\
        .to_csv(tablename, sep='\t',index=False)
    print 'Written: ', tablename

    # write bedtools intersect table
    out_bed_name = ref_path + '/methyl_table.bed'
    coor_names = ['chrom','start','end','Genomic location','Tissue']
    df = df[coor_names]
    df.to_csv(out_bed_name,sep='\t', 
            header=False, index=False)
    print 'Written: ', out_bed_name


if __name__ == '__main__':
    main()

