#!/bin/bash

BAMPATH=$1
SAMPLENAME=$2
RESULTPATH=$3
MEMORY=$4
PICARD=/work/02727/cdw2854/lonestar/src/broadinstitute-picard-b7a1335/dist/picard.jar

if [ $# != 4 ]
then
	echo usage: $0 \<bampath\> \<samplename\> \<resultpath\> \<memory\>
	exit
fi


FILES=$(ls $BAMPATH/$SAMPLENAME* | tr '\n' ' ' | rev |cut -d' ' -f2- | rev )
FILES=${FILES// / I=}
java -Xmx$MEMORY -jar $PICARD GatherBamFiles I=$FILES \
				O=/dev/stdout \
| java -Xmx$MEMORY -jar $PICARD  SortSam \
				INPUT=/dev/stdin \
				OUTPUT=/dev/stdout 	\
				SORT_ORDER=coordinate \
| java -Xmx$MEMORY -jar $PICARD MarkDuplicatesWithMateCigar \
			ASSUME_SORTED=true \
			REMOVE_DUPLICATES=true \
			INPUT=/dev/stdin \
			METRICS_FILE=$RESULTPATH/$SAMPLENAME.txt \
			OUTPUT=/dev/stdout \
> $RESULTPATH/$SAMPLENAME.bam 
echo Finished merging $SAMPLENAME
