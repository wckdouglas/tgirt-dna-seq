#!/bin/env python

import pandas as pd
import os
import MySQLdb as db
import numpy as np

def sqlTotable(tableResult, cleantable, sourceDict):
    for record in tableResult:
        chrom, end, start = record['chrom'], int(record['chromEnd']), int(record['chromStart'])
        coordinate = '%s\t%i\t%i' %(chrom.replace('chr',''), start, end)
        scores = record['sourceScores'].split(',')[:-1]
        ids = record['sourceIds'].split(',')[:-1]
        assert len(ids) == len(scores), 'Wrong id and score numbers'
        [cleantable.write(coordinate + '\t' + sourceDict[sourceID] + '\t' + sourceScore + '\n') \
            for sourceID, sourceScore in zip(ids,scores)]
    return 0

def UCSCsqlDB(ucsc):
    ucsc.query('SELECT * FROM wgEncodeRegDnaseClustered;')
    result = ucsc.store_result().fetch_row(how=1, maxrows=0)
    return result

def makingIDdict(ucsc):
    ucsc.query('SELECT * FROM wgEncodeRegDnaseClusteredSources;')
    result = ucsc.store_result().fetch_row(how=1, maxrows=0)
    sourceDict = {str(record['id']): record['name'] for record in result}
    return sourceDict

def cleaningBed(ucsc,sourceDict, cleanTablename):
    tableResult = UCSCsqlDB(ucsc)
    with open(cleanTablename,'w') as cleantable:
        cleantable.write('chrom\tstart\tend\tsource\tscore\n')
        sqlTotable(tableResult, cleantable, sourceDict)
    return cleanTablename


def readTable(cleanTablename, bedFilename, cellTablename):
    df = pd.read_csv(cleanTablename,sep='\t')
    cellDF = pd.read_csv(cellTablename, sep='\t',header=None).iloc[:,[1,3,4,5]]
    cellDF.columns = ['source','tissue','cellType','type']
    df = pd.merge(df,cellDF,how='inner',on='source')
    df = df[(df['type']=='normal')]# & (df['cellType']=='primary cell')]
    df.drop(['cellType','type','source'],axis=1,inplace=True)
    df = df.groupby(['chrom','start','end','tissue']).sum().reset_index()
    df = pd.pivot_table(df,values='score',columns='tissue',index=['chrom','start','end'],aggfunc=np.mean)
    df.reset_index(inplace=True)
    df.fillna(0,inplace=True)
    annotation = df.iloc[:,:3]
    data = df.iloc[:,3:]
    #data = data.apply(lambda x: np.true_divide(x,max(x)))
    annotation['name'] = ['cluster_%i' %i for i in np.arange(len(annotation))]
    annotation['score'] = np.repeat('0', len(annotation))
    annotation['strand'] = np.repeat('+',len(annotation))
    df = pd.concat([annotation, data], axis=1)
    df.to_csv(bedFilename,sep='\t', index=False)
    print 'Written %s' %bedFilename

def main():
    datapath = '/scratch/02727/cdw2854/plasmaDNA/ctcfData'
    datapath = '/Users/wckdouglas/plasmaDNA/reference/ctcfData'
    datapath = '/scratch/cdw2854/plasmaDNA/CTCFdata'
    cleanTablename = datapath + '/cleanedCTCF.tsv'
    bedFilename = datapath + '/cellCTCF.bed'
    cellTablename = datapath + '/cellType.tsv'
    if not os.path.isfile(cellTablename):
        os.system('curl http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/wgEncodeCell.txt.gz |'+
             'gunzip > %s' %cellTablename)
        print 'downloaded cell type description'
    if not os.path.isfile(cleanTablename):
        ucsc = db.connect(host='genome-mysql.cse.ucsc.edu',user='genome',db="hg38")
        sourceDict = makingIDdict(ucsc)
        print 'Made index'
        cleanTablename = cleaningBed(ucsc, sourceDict, cleanTablename)
        ucsc.close()
        print 'Written %s' %cleanTablename
    print 'Using %s' %cleanTablename
    readTable(cleanTablename, bedFilename, cellTablename)

if __name__ == '__main__':
    main()

