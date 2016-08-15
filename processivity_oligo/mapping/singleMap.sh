#!?bin/bash

PROJECTPATH=/scratch/02727/cdw2854/TGIRT_plasma_dna/syntheticOligo
DATAPATH=$PROJECTPATH/errorFreeFastq
RESULTPATH=$PROJECTPATH/bamFile
INDEX=/corral-repl/utexas/2013lambowitz/Ref/syntheticOligos/celSpikeinLong/spikein.fa
THREADS=12

for FQ in $DATAPATH/*fq.gz
do
    SAMPLENAME=$(basename ${FQ%.fq.gz})
    echo bwa mem -t $THREADS \
        $INDEX \
        $FQ \
    \| samtools view -b@ $THREADS - \
    \| samtools sort \
        -@ $THREADS  \
        -T $RESULTPATH/$SAMPLENAME \
        -O bam - \
    \> $RESULTPATH/$SAMPLENAME.bam 
done
