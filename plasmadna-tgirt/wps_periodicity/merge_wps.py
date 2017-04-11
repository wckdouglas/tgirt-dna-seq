#!/bin/bash

import pandas as pd
import glob

work_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS'
analytics = ['periodicity_gene_body', 'periodicity_tss']
ffts = ['R','scipy']
for analytic in analytics:
    for fft in ffts:
        files = glob.glob(work_path + '/' + analytic + '/' + fft + '/*.bed')
        samplenames = set(map(lambda x: x.split('.')[0], files))
        for samplename in list(samplenames):
            out_table_name = samplenames + '.tsv'
            table_names = glob.glob(samplename+'*.bed')
            df = map(read_table, table_names)
            df = pd.concat(df, axis=0)
            df.to_csv(out_table_name, sep='\t', index=False)
            print 'Written: ', out_table_name
