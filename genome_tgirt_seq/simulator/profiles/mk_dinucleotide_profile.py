#!/usr/bin/env python

import sys
import pandas as pd
import os


if len(sys.argv) != 2:
    sys.exit('usage: python %s <tablename>' %sys.argv[0])
tablename = sys.argv[1]
out_table_name = os.path.basename(tablename).split('.')[0] + '.csv'


pd.read_table(tablename) \
        .pipe(pd.melt, id_vars = ['total','index','End','sample'],
                var_name = 'dinucleotide', value_name = 'fraction') \
        .pipe(lambda d: d[~d.dinucleotide.str.contains('N')]) \
        .query('index==-1') \
        .drop(['sample','total','index'],axis=1)\
        .to_csv(out_table_name,index=False)
