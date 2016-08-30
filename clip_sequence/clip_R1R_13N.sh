#!/bin/bash

DATAPATH=$DATA/JA16493
CLIPPED_PATH=$DATAPATH/clipped
mkdir -p $CLIPPED_PATH

for FQ1 in $DATAPATH/[g]*R1_001.fastq.gz
do
	R1_NAME=$(basename $FQ1)
	FQ2=${FQ1/R1/R2}
	R2_NAME=$(basename $FQ2)

	#Read 1 process
	echo seqtk trimfq -b18 $FQ1 \| gzip \> $CLIPPED_PATH/$R1_NAME
	if [[ $R1_NAME == DB* ]]
	then
		echo seqtk trimfq -b41 $FQ2 \| gzip \> $CLIPPED_PATH/$R2_NAME
	else
		echo cp $FQ2 $CLIPPED_PATH/$R2_NAME
	fi
done
