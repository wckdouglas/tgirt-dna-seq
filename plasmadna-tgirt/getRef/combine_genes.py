#!/usr/bin/env python

import pandas as pd
import glob 

file_path = '/stor/work/Lambowitz/cdw2854/plasmaDNA/genes/'
files = glob.glob(file_path+'*.csv.zip')
df = map(lambda x: pd.read_csv(x), files)
df = pd.concat(df,axis=0)
out_file = file_path + '/rna.csv'
df.to_csv(out_file, index=False)
print 'Made %s' %out_file

