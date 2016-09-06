#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasma_project
DATAPATH=$PROJECTPATH/rawData
INDEXPATH=/corral-repl/utexas/2013lambowitz/Ref/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA.fa
ADAPTORS=adaptors.fa
CORES=12
MEMORY=1g

for FQ in $DATAPATH/PDB*R1_001.fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	echo python DNAmapping.py --fq1=$FQ --outdir=$PROJECTPATH --index=$INDEX \
						--threads=$CORES --adaptor=$ADAPTORS \
						--memory=$MEMORY
done 
