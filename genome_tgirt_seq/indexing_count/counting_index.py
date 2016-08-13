#!/usr/bin/env python

import pandas  as pd
import cjson
import numpy as np
import glob
from functools import partial
from multiprocessing import Pool
import matplotlib
matplotlib.use('Agg')
import seaborn as sns
import matplotlib.pyplot as plt
import re
sns.set_style('white')

def extractJSON(json_line):
    record = cjson.decode(json_line)
    barcode = record[0]
    table = record[1]
    seq_count = len(table)
    return barcode, seq_count

def writeLine(json_file, table_name):
    with open(json_file,'r') as json_handle, open(table_name,'w') as out_file:
        if 'double' in json_file:
            out_file.write('barcode\tdouble_count\tindex1\tindex2\n')
            for json_line in json_handle:
                barcode, seq_count = extractJSON(json_line)
                index1, index2 = barcode.split('|')
                out_file.write('%s\t%i\t%s\t%s\n' %(barcode, seq_count, index1, index2))
        else:
            out_file.write('index1\tsingle_count\n')
            for json_line in json_handle:
                barcode, seq_count = extractJSON(json_line)
                out_file.write('%s\t%i\n' %(barcode, seq_count))
    return 0

def analyzeJSON(result_path,json_file):
    samplename = json_file.split('/')[-1].split('.')[0]
    print 'Analyzing %s' %samplename
    table_name = result_path + '/' + samplename + '.tsv'
    writeLine(json_file,table_name)
    print 'Written %s' %table_name
    return table_name

def compareFile(double_file):
    single_file  = double_file.replace('-double-BC','')
    df_d = pd.read_table(double_file)\
        .drop('barcode',axis=1)

    df = pd.read_table(single_file)\
        .merge(df_d, how='inner',on='index1') \
        .groupby(['single_count','double_count'])\
        .size() \
        .rename('count') \
        .reset_index()
    return df

def combineTables(table_names):
    # combining tables
    double_bc_table_names = filter(lambda x: re.search('double-BC',x), table_names)
    p = Pool(12)
    df = p.map(compareFile,double_bc_table_names)
    p.close()
    p.join()
    df = pd.concat(df,axis=0) \
        .groupby(['single_count','double_count'])\
        .agg({'count':'sum'}) \
        .reset_index()
    return df

def plotDensity(df, figurename):
    with sns.plotting_context('paper',font_scale=1.3):
        p = sns.FacetGrid(data = df)
    p.map(sns.kdeplot, 'single_count','double_count',weights='density')
    p.add_legend()
    p.set_xlabels('Single barcode family')
    p.set_ylabels('Double barcode family')
    p.savefig(figurename)

def main():
    json_path = '/stor/work/Lambowitz/Data/NGS/JA16381/combined/splitted'
    result_path = '/stor/work/Lambowitz/cdw2854/jurkatCells/index_tables'
    figurename = result_path + '/barcode_count.pdf'
    json_files = glob.glob(json_path + '/DB*.json')
    func = partial(analyzeJSON, result_path)

    # make tables
    p = Pool(12)
    table_names = p.map(func, json_files)
    p.close()
    p.join()

    df = combineTables(table_names)
    df.to_csv(result_path + '/barcode_count.tsv',sep='\t',index=False)
    plotDensity(df,figurename)
    return 0

if __name__ == '__main__':
    main()
