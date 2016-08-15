#!/bin/bash

PROJECTPATH=/stor/work/Lambowitz/cdw2854/processivity
FASTQPATH=$PROJECTPATH/splitted_data
BAMPATH=$PROJECTPATH/bamFile
INDEX=$REF/syntheticOligos/celSpikeinLong/spikein.fa
THREADS=6
SUFFIX=R1_001.fastq.gz
mkdir -p $BAMPATH

for FQ1 in $FASTQPATH/*$SUFFIX
do
    SAMPLE=$(basename ${FQ1%_$SUFFIX})
    FQ2=${FQ1/R1/R2}
    echo bwa mem -t $THREADS $INDEX $FQ1 $FQ2 \
		\| samtools view -b@ $THREADS \
		\| bamtools filter -script flag_filter.json \
        \| samtools sort -@ $THREADS -O bam -T $BAMPATH/$SAMPLE - \
        \> $BAMPATH/${SAMPLE}.bam 
done
