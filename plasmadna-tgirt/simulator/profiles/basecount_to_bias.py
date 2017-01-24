import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

df = pd.read_table('refD3_00_S10_1P.head.baseCount.tsv')
df['rowsum']=df[list('ACTG')].values.sum(axis=1)


for base in list('ACTG'):
    df[base] = np.true_divide(df[base],df['rowsum'])
df['index'] = np.arange(len(df))
df.drop('rowsum',axis=1)
df.to_csv('bias.csv',index=False)

