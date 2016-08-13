#!/usr/bin/env python

from pybedtools import BedTool
import seaborn as sns
import pandas as pd
import numpy as np
from collections import defaultdict
import os
sns.set_style('white')

def percentage(x):
    return np.true_divide(x,np.sum(x))*100

def plot_coverage(df, figurename):
    with sns.plotting_context('paper',font_scale=1.4):
	p = sns.FacetGrid(data = df)
    p.map(sns.barplot, 'coverage','position_count')
    p.set(xlabel='Coverage', ylabel='Percentage of positions')
    p.savefig(figurename)
    print 'Saved: %s' %figurename
    return 0

def parse_bam_depth(bam_file, genome_file):
    samplename = os.path.basename(bam_file).split('.')[0]
    position_cov = defaultdict(int)
    print 'Running %s' %samplename
    for position in BedTool().genome_coverage(d=True, g=genome_file, ibam=bam_file):
	fields = position.fields
	coverage = int(fields[2])
	position_cov[coverage] += 1
    print 'Finished %s' %samplename
    df = pd.DataFrame({'coverage': position_cov.keys(),
		    'position_count':position_cov.values()})\
	.assign(position_count_percentage = lambda x: percentage(x['position_count']))
    return df

def main():
    project_path = '/scratch/02727/cdw2854/jurkatCells'
    ref_path = '/corral-repl/utexas/2013lambowitz/Ref/hg19/Sequence/WholeGenomeFasta'
    figure_path = project_path + '/figures'
    genome_file = ref_path + '/hg19.genome'
    bam_path = project_path + '/merged_bam'
    bam_file = bam_path + '/jurkat_DB_SB.sorted.bam'
    figurename = figure_path + '/coverage.pdf'
    table_name = bam_path + 'coverage.tsv'
    df = parse_bam_depth(bam_file, genome_file)
    df.to_csv(table_name,index=False,sep='\t')
    plot_coverage(df,figurename)
    return 0


if __name__ == '__main__':
    main()
