#!/usr/bin/env python

import matplotlib
matplotlib.use('Agg')
from pybedtools import BedTool
import glob
import os
import pandas as pd
import numpy as np
import pylab as plt
from multiprocessing import Pool
import seaborn as sns
import dmisc.folder_action as folders
sns.set_style('white')

ncRNA = ["sense_intronic","3prime_overlapping_ncrna",'processed_transcript',
        'sense_overlapping','Other_lncRNA']
smncRNA = ['misc_RNA','snRNA','piRNA']
large_rRNA = ['28S_rRNA','18S_rRNA']
small_rRNA = ['rRNA','5S_rRNA','58S_rRNA','5.8S_rRNA']
protein_coding = ['protein_coding','TR','IG']
def groupGeneType(x):
    rnatype = 'Other ncRNA' if x in ncRNA else \
                  'TR'  if 'TR' in x else \
                  'IG' if 'IG' in x else \
                  'Mt' if 'Mt_' in x else \
                  'tRNA' if 'tRNA' in x else \
                  '5/5.8S rRNA' if x in small_rRNA else \
                  '18/28S rRNA' if x in large_rRNA else \
                  'Other sncRNA'if x in smncRNA else \
                  'Pseudogenes' if 'pseudogene' in x  else x
    return rnatype

def intersect2dataframe(intersection, samplename):
    scores = []
    rnaTypes = []
    ids = []
    names = []
    for feature in intersection:
        fields = feature.fields
        score = float(fields[6]) 
        rnaType = fields[13]
        name = fields[10]
        gene_id = fields[-1]
        ids.append(gene_id)
        scores.append(score)
        rnaTypes.append(rnaType)
        names.append(name)
    df = pd.DataFrame({'id':ids, 'name': names, 'scores': scores, 'rnaType': rnaTypes})
    df['samplename'] = np.repeat(samplename, len(df))
    return df

def findRNA(args):
    bedfile, rnaBed, motifBed = args
    samplename = os.path.basename(bedfile).split('.')[0]
    print 'Running %s' %samplename
    shortBed = BedTool(bedfile)
    longBed = BedTool(bedfile.replace('Short','Long'))
    intersection = shortBed.intersect(b=motifBed, f=0.5, r = True, v=True) \
            .intersect(b=longBed, v=True)\
            .intersect(b=rnaBed, f = 0.5, wb=True, s=True) 
    df = intersect2dataframe(intersection, samplename)
    return df

def returnPercentage(df):
    df['scores'] = np.true_divide(df['scores'],np.sum(df['scores'])) *100
    return df

def plotType(df, figurename):
    df['rnaType'] = map(groupGeneType, df['rnaType'])
    df = df.groupby(['rnaType','samplename']).sum().reset_index()
    df = df.groupby('samplename').apply(returnPercentage)
    df = df.pivot_table(values='scores', columns = 'rnaType', index='samplename', aggfunc=np.sum)
    fig = plt.figure(figsize=(18, 18))
    p = df.plot.bar(stacked=True, colormap='Pastel2_r')
    lgd = p.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    p.set_ylim(0,100)
    plt.savefig(figurename, dpi=300, bbox_extra_artists=(lgd,), bbox_inches='tight')
    print 'Saved %s' %figurename
    return 0

def main():
    projectpath = '/stor/work/Lambowitz/cdw2854/plasmaDNA'
    bedpath = projectpath + '/genomeWPS/unstranded'
    bedpath = projectpath + '/genomeWPS/stranded'
    resultpath = projectpath + '/rnaCalls'
    referencepath = '/stor/work/Lambowitz/cdw2854/reference/hg38'
    bedfiles = glob.glob(bedpath + '/*bed')
    rnaRef = referencepath + '/genes_count.bed'
    motifRef = referencepath + '/MotifFeatures.gff'
    resultFilename = resultpath + '/rnaCount.tsv'
    figurename = resultpath + '/rnaCount.pdf'
    folders.makeFolder(resultpath)

    #Run files 
    rnaBed = BedTool(rnaRef)
    motifBed = BedTool(motifRef)
    pool = Pool(24)
    dfs = pool.map(findRNA, [(bedfile, rnaBed, motifBed) for bedfile in bedfiles])
    df = pd.concat(dfs, axis=0)
    df.to_csv(resultFilename, sep='\t', index=False)
    print 'Written %s' %resultFilename
    p = plotType(df, figurename)
    plotDF(df, figurename)
    return 0

if __name__ == '__main__':
    main()
