#!/bin/bash

DATAPATH=$DATA/JA16466
CLIPPED_PATH=$DATA/clipped
mkdir -p $CLIPPED_PATH

for FQ1 in $DATAPATH/[PD]*R1_001.fastq.gz
do
	R1_NAME=$(basename $FQ1)
	FQ2=${FQ1/R1/R2}
	R2_NAME=$(basename $FQ2)

	#Read 1 process
	echo zcat $FQ1 \| fastx_trimmer -f 19 \| gzip \> $CLIPPED_PATH/$R1_NAME
	if [[ $R1_NAME == DB* ]]
	then
		echo zcat $FQ2 \| fastx_trimmer -f 42 \| gzip \> $CLIPPED_PATH/$R2_NAME
	else
		echo cp $FQ2 $CLIPPED_PATH/$R2_NAME
	fi
done
