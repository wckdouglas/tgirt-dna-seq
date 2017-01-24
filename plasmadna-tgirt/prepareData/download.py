#!/usr/bin/env python

import pandas as pd
import os

df = pd.read_table('SraRunTable.txt', sep='\t') \
    .pipe(lambda d: d[d.disease_s.str.contains('Healthy|Breast')])

out_dir = os.environ['WORK'] + '/cdw2854/plasmaDNA/raw_data'
base_url = 'ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR'
for run in df.Run_s:

    if '993' not in run:
        srr_prefix = run[:6]
        command = 'curl -o {out_dir}/{run}.sra {base}/{prefix}/{run}/{run}.sra'.format(base=base_url,
                                                  prefix = srr_prefix,
                                                  out_dir = out_dir,
                                                  run = run)
        print command



