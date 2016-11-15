#!/bin/bash

TARGET_DIR=/stor/work/Lambowitz/cdw2854/target-seq/raw_reads
FROM="$DATA/JA16637 $DATA/SA16172/JA16594"
for DIR in $FROM
do
	for FQ in $DIR/*fastq.gz
	do
		FILE_NAME=$(basename $FQ)
		echo ln -s $FQ $TARGET_DIR/$FILE_NAME
	done
done
