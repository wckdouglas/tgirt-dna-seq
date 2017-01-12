#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
DATAPATH=$PROJECTPATH/rawData/umi2id
RESULTPATH=$PROJECTPATH/splittedFastq
SRCPATH=/work/02727/cdw2854/lonestar/src
PROGRAM=$SRCPATH/fastq-tools/bin/splitFastq
mkdir -p $RESULTPATH

for FASTQ in $(ls $DATAPATH/*gz | grep -v 'try\|93\|50')
do 
	SAMPLE=$(basename $FASTQ)
	SUFFIX=$(echo $SAMPLE | rev | cut -d'.' -f1 | rev )
	if [[ $SUFFIX == gz ]]
	then
		SUFFIX=$(echo $SAMPLE | rev | cut -d'.' -f1-2 | rev )
	fi
	SAMPLENAME=${SAMPLE%.$SUFFIX}
    echo time $PROGRAM -i $FASTQ -n 5000000 \
		-o $RESULTPATH/$SAMPLENAME -z 
done
