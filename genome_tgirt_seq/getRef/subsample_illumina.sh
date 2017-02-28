#!/usr/bin/bash

PROJECTPATH=/stor/work/Lambowitz/cdw2854/genomeDNA
IN_PATH=$PROJECTPATH/rawData
OUT_PATH=$PROJECTPATH/bam_files/clustered_fastq

for i in $(seq 1 3)
do
	echo sambamba view  \
		--subsample=0.006 \
		--subsampling-seed=${i} \
		--format=bam \
		$IN_PATH/NA12878_S1.nameSorted.bam \
	\| bedtools bamtofastq -i - \
		-fq /dev/stdout -fq2 /dev/stdout \
	\| seqtk dropse  - \
	\| gzip \
	\> $OUT_PATH/NA12878_illumina_${i}.fastq.gz
done

