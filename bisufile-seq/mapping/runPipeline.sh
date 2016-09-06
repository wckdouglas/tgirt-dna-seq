#!/bin/bash

PROJECTPATH=/stor/work/Lambowitz/cdw2854/bisufite_seq
DATAPATH=$DATA/JA16234/bisulfite
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA.fa
ADAPTORS=adaptors.fa
CORES=24
MEMORY=1g

for FQ in $DATAPATH/PDB*R1_001.fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	echo python DNAmapping.py --fq1=$FQ --outdir=$PROJECTPATH --index=$INDEX \
		--threads=$CORES --adaptor=$ADAPTORS \
		--memory=$MEMORY
done 
