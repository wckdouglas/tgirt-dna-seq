#!/bin/bash

PROJECTPATH=/stor/work/Lambowitz/cdw2854/bisufite_seq/
DATAPATH=$PROJECTPATH/raw_Data
DATAPATH=$DATAPATH/UMI_clustered/clustered_fastq
INDEXPATH=$REF/hg19/Sequence/WholeGenomeFasta
INDEX=$INDEXPATH/genome.fa
ADAPTORS=adaptors.fa
CORES=24

for FQ in $DATAPATH/P13**R1_001.fastq.gz
do
	SAMPLENAME=$(basename $FQ)
	if [[ $SAMPLENAME == *B* ]]
	then
		echo python DNA_meth_mapping.py --fq1=$FQ \
			--outdir=$PROJECTPATH --index=$INDEX \
			--threads=$CORES --adaptor=$ADAPTORS 
	else
		echo python DNA_mapping.py --fq1=$FQ \
			--outdir=$PROJECTPATH --index=$INDEX \
			--threads=$CORES --adaptor=$ADAPTORS 
	fi
done 
