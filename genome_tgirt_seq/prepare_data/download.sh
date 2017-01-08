TARGET_DIR=/stor/work/Lambowitz/cdw2854/ecoli_genome/rawData/k12
#LINK=" ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R1.fastq.gz  \
#		ftp://webdata:webdata@ussd-ftp.illumina.com/Data/SequencingRuns/MG1655/MiSeq_Ecoli_MG1655_110721_PF_R2.fastq.gz \
#		http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/2015.07.07/Ecoli_S1_L001_R1_001.fastq.gz \
#		http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/2015.07.07/Ecoli_S1_L001_R2_001.fastq.gz "
LINK=http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/2015.07.07/Ecoli_S1_L001_R2_001.fastq.gz 
for FILE in $LINK
do
	SAMPLENAME=$(echo $FILE | rev | cut -d'/' -f1 | rev)
	echo curl  $FILE \| seqtk trimfq -L 150 - \| gzip \> $TARGET_DIR/$SAMPLENAME
done

SRR='SRR522163 SRR519925 SRR519926'
for SR in $SRR
do
	echo fastq-dump --split-3 --gzip -O $TARGET_DIR $SR
done

for SRR in $(cat files.txt)
do
	echo fastq-dump --split-3 --gzip -O $TARGET_DIR $SRR
done

