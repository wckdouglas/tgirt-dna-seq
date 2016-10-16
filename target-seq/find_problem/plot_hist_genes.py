#/usr/bin/env python

import seaborn as sns
import matplotlib.pyplot as plt
import pandas as pd
from pybedtools import BedTool, Interval, set_tempdir
import glob
import numpy as np
from functools import partial
from multiprocessing import Pool
import pysam
import re
sns.set_style('white')

def bedpeToBed(bedline):
    fields = bedline.fields
    start = np.min(map(long,[fields[1],fields[4]]))
    end = np.max(map(long,[fields[2],fields[5]]))
    chrom = fields[0]
    name = fields[6]
    strand = fields[-2] 
    length = end-start
    return Interval(chrom=chrom,
                    start = start,
                    end = end,
                    strand = strand,
                    score=length,
                    name = name)

def filterBam(in_bam):
    out_bam = in_bam.replace('.bam','.filtered.bam')
    out_count = 0
    with pysam.Samfile(in_bam,'rb') as inbam:
        with pysam.Samfile(out_bam,'wb',template=inbam) as outbam:
            for aln in inbam:
                if not aln.is_unmapped and aln.isize != 0:
                    outbam.write(aln)
                    out_count += 1
    print 'Exported %i alignments' %(out_count/2)
    return out_bam

def bamToBed(bam_file):
    bed_file_name = bam_file.replace('.bam','.bed')
    bed = BedTool(bam_file) \
        .bam_to_bed(bedpe=True, mate1=True) \
        .filter(lambda f: f.fields[0]==f.fields[3])\
        .filter(lambda f: f.fields[-1]!=f.fields[-2])\
        .each(bedpeToBed)\
        .filter(lambda f: f.score > 200)\
        .sort() \
        .moveto(bed_file_name)
    return bed

def coverageDF(gene_df, hist_bed, bam_file, samplename):
    bed = bamToBed(bam_file)
    cov_df = BedTool(hist_bed)\
        .coverage(b=bed,d=True)\
        .to_dataframe(names = ['id','start','end','position','depth']) \
        .assign(sample_id = samplename) \
        .merge(gene_df)
    return cov_df

def processBam(hist_bed, gene_df, bam_file):
    samplename = bam_file.split('/')[-1].split('.')[0]
    print 'Running %s '%samplename
    filtered_bam = filterBam(bam_file)
    print 'Filtered %s' %filtered_bam
    cov_df = coverageDF(gene_df, hist_bed, filtered_bam, samplename)
    figurename = bam_file.replace('.bam','.pdf')
    with sns.plotting_context('paper',font_scale = 1.2):
        p = sns.FacetGrid(data = cov_df, col = 'name', 
                col_wrap=4, sharex=False, sharey=False)
    p.map(sns.barplot, 'position', 'depth',color='steelblue')
    p.set_xticklabels(rotation = 70) 
    p.savefig(figurename)
    print 'Saved: %s' %figurename

def main():
    project_path = '/stor/work/Lambowitz/cdw2854/target-seq'
    tx_path = '/stor/work/Lambowitz/ref/human_transcriptome'
    bam_path = project_path + '/name_sorted'
    tx_len = tx_path + '/transcript_count.bed'
    tx_info = tx_path + '/tx.info'
    gene_info = '/stor/work/Lambowitz/ref/GRCh38/Bed_for_counts_only/genes.Info'
    set_tempdir(bam_path)

    bam_files = glob.glob(bam_path + '/*bam')
    bam_files = filter(lambda x: not re.search('filtered',x), bam_files)

    gene_df = pd.read_table(gene_info) \
        .rename(columns = {'GeneID':'gene'}) \
        .drop(['type'],axis = 1)\
        .merge(pd.read_table(tx_info), on='gene',how='inner') \
        .pipe(lambda d: d[d['name'].str.contains('HIST1H3')])

    hist1h3_genes = pd.read_table(tx_len, names = ['id','start','end'])\
        .pipe(lambda d: d[np.in1d(d['id'],gene_df['id'])])

    hist_bed = BedTool()\
        .from_dataframe(hist1h3_genes)

    bamFunc = partial(processBam, hist_bed, gene_df)

    Pool(4).map(bamFunc,bam_files)

if __name__ == '__main__':
    main()
