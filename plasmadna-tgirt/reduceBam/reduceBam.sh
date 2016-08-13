#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BAMPATH=$PROJECTPATH/bamFiles
MERGEBAM=$PROJECTPATH/mergedBam
MEMORY=2g

if [ -d $PROJECTPATH ]
then
	mkdir -p $MERGEBAM
else
	echo $PROJECTPATH not exists
	exit
fi


for SAMPLE in $(ls $BAMPATH/*bam | rev |cut -d'/' -f1 | rev |  cut -d'_' -f1 | sort | uniq )
do
	echo sh mergeBam.sh $BAMPATH $SAMPLE $MERGEBAM $MEMORY
done
