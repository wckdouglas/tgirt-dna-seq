#!/bin/bash

PROJECTPATH=$WORK/cdw2854/ecoli_genome
#PROJECTPATH=$WORK/cdw2854/genomeDNA
DATAPATH1=$PROJECTPATH/rawData/k12
DATAPATH2=$PROJECTPATH/rawData/k12/CUT_75BP
INDEXPATH=$REF/Ecoli
INDEX=$INDEXPATH/k12_mg1655.fa
#INDEX=$REF//hg19/Sequence/WholeGenomeFasta/genome.fa
PYTHON=$(which python)
ADAPTORS=adaptors.fa
CORES=1

for DATAPATH in $DATAPATH1 $DATAPATH2
do
	for FQ in $DATAPATH/*R1_001.fastq.gz 
	do
		SAMPLENAME=$(basename $FQ)
		if [[ $SAMPLENAME == K12* ]]
		then
			ADAPTORS=adaptors.fa
		else
			ADAPTORS=TruSeq2-PE.fa
		fi
		echo $PYTHON DNAmapping.py --fq1=$FQ  \
			--outdir=$PROJECTPATH --index=$INDEX \
			--threads=$CORES --adaptor=$ADAPTORS 
	done
done 
