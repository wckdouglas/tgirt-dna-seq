#!/bin/bash

TARGET_DIR=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12
FROM="$DATA/JA16703 $DATA/JA16720 $DATA/JA16791 $DATA/JA17030" 
for DIR in $FROM
do
	for FQ in $DIR/K12*fastq.gz
	do
		FILE_NAME=$(basename $FQ)
		ln -s $FQ $TARGET_DIR/$FILE_NAME
		echo linked $FQ to $TARGET_DIR/$FILE_NAME
	done
done
