#!/bin/bash

TARGET_DIR=/stor/work/Lambowitz/cdw2854/target-seq/raw_reads
FROM="$DATA/JA16637 $DATA/SA16172/JA16594 $DATA/JA16663 $DATA/JA16810"
for DIR in $FROM
do
	for FQ in $DIR/*fastq.gz
	do
		FILE_NAME=$(basename $FQ)
		ln -s $FQ $TARGET_DIR/$FILE_NAME
	done
done
