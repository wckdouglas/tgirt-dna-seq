#!/bin/bash

DESTINATION=/scratch/02727/cdw2854/TGIRT_plasma_dna/plasmaResults/data
DBPATH=/corral-repl/utexas/Exosome-RNA-seq/NCBI|
DBLINK=http://dl.dropbox.com/u/51653511/SRAmetadb.sqlite.gz
CURRENT=`pwd`

#cd $DBPATH 
#wget $DBLINK
#gunzip SRAmetadb.sqlite.gz

cd $DESTINATION
Rscript $CURRENT/download.R
for LINK in `cat filelist.tsv | egrep 'SRR2130050|SRR2130051|SRR2129993|SRR2130052' | cut -f5`
do
   curl -O $LINK 
done
