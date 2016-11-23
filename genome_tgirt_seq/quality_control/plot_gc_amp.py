import numpy as np
import pandas as pd

def autocorr(x): 
    result = np.correlate(x, x, mode='full') 
    return result[result.size/2:]

df = pd.read_table('CTCFwps.tsv') 

wps = df\
    .pipe(lambda d: d[d['type'].str.contains('Long')])\
    .pipe(lambda d: d[d['samplename']=='SRR2130051'])\
    .wps
    
ax = plt.subplot(111)
plt.plot(autocorr(wps))
plt.xlim(0,300)
ax.set_xticks(range(0,301,10))
ax.set_xticklabels(range(0,301,10), rotation=90)