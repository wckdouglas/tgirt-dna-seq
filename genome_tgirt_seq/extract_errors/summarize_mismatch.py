#!/usr/bin/env python

import pyximport
pyximport.install()
from parse_pilebed import summarize_mismatch
import glob
from multiprocessing import Pool
import pandas as pd
import sys

mismatch_table_path = '/stor/work/Lambowitz/cdw2854/ecoli_genome/mismatch_profiles'
files = glob.glob(mismatch_table_path + '/*tsv')
p = Pool(24)
df = p.map(summarize_mismatch ,files)
p.close()
p.join()
df = pd.concat(df)
result_table = mismatch_table_path + '/mismatch_profiles.tsv'
df.to_csv(result_table, sep='\t', index=False)
sys.stderr.write('Written %s\n' %(result_table))
