#!/bin/bash

PROJECTPATH=/scratch/02727/cdw2854/plasmaDNA
BAMPATH=$PROJECTPATH/bamFiles
BEDPATH=$PROJECTPATH/dupbedFiles
BEDTOBEDPE=$(which bedpeTobed)
BEDTOOLS=/work/02727/cdw2854/lonestar/src/bedtools2/bin/bedtools
PYTHON=$(which python)
BAMTOOLS=$(which bamtools)
SAMTOOLS=$(which samtools)
PICARD_TOOLS=/work/02727/cdw2854/lonestar/src/broadinstitute-picard-b7a1335/dist/picard.jar
mkdir -p $BEDPATH

for SAMPLENAME in `ls $BAMPATH/*bam | rev | cut -d'/' -f1 | rev | cut -d'_' -f1 | sort | uniq`
do
	TEMP_DIR=$BEDPATH/${SAMPLENAME}_tmp
	TMP_FILE=$TEMP_DIR/$SAMPLENAME
	mkdir -p $TEMP_DIR
	echo $SAMTOOLS cat $BAMPATH/$SAMPLENAME*bam \
		\| $BAMTOOLS filter -in - -script filterCriteria.json  \
		\| $SAMTOOLS sort -T ${TMP_FILE}.coordincate -O bam  \
		\| $SAMTOOLS rmdup - - \
		\| $SAMTOOLS sort -T ${TMP_FILE}.name -O bam -n - \
		\| $BEDTOOLS bamtobed -i - -mate1 -bedpe \
		\| $BEDTOBEDPE -i - \
		\> $BEDPATH/${SAMPLENAME}.bed  
done
