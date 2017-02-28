BAM=ftp://ftp.sra.ebi.ac.uk/vol1/ERA172/ERA172924/bam/NA12878_S1.bam
DESTINATION=/stor/work/Lambowitz/cdw2854/genomeDNA/rawData

curl $BAM \
	> $DESTINATION/NA12878_S1.bam
#	| sambamba sort  \
#		--tmpdir=$DESTINATION/temp \
#		--nthreads 12 \
#		--memory-limit 10GB \
#		 --sort-by-name \
#		 /dev/stdin \
