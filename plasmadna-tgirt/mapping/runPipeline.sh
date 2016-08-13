#!/bin/bash

PROJECTPATH=$SCRATCH/plasmaDNA
DATAPATH=$PROJECTPATH/rawData/clipped
INDEXPATH=$REF/GRCh38/hg38_rDNA
INDEX=$INDEXPATH/genome_rDNA
PYTHON=$(which python)
CORES=12

for FQ in `ls $DATAPATH/*R1_001.fastq.gz $DATAPATH/*R1_001.fq.gz | grep 'GJ\|RNA\|ref' -v`
do
	SAMPLENAME=$(basename $FQ)
	if [[ $SAMPLENAME == SRR* ]]
	then
		ADAPTORS=TruSeq2-PE.fa
	#elif [[ $SAMPLENAME == DB* ]]
	#then
	#	ADAPTORS=double_indexed_adaptors.fa
	#elif [[ $SAMPLENAME == *errorFree* ]]
	#then
	#	ADAPTORS=indexed_adaptors.fa
	else
		ADAPTORS=adaptors.fa
	fi
	echo $PYTHON DNAmapping.py --fq1=$FQ --outdir=$PROJECTPATH --index=$INDEX \
						--threads=$CORES --adaptor=$ADAPTORS 
done 
