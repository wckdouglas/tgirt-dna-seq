#!/usr/bin/env python

import pandas as pd
import seaborn as sns
import glob

filepath = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genomeWPS/periodicity_tables'
filenames = glob.glob(filepath + '/*tsv')
dfs = map(pd.read_table, filenames)
df = pd.concat(dfs, axis = 0)  \
    .drop(['chrom','start','end'], axis = 1) \
    .pivot(index = 'name', values = 'periodicity', columns = 'samplename') \
    .reset_index() 

np.corrcoef(df.drop(['name'],axis = 1))



