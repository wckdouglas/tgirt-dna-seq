#!/bin/bash

FASTQPATH=/scratch/02727/cdw2854/TGIRT_plasma_dna/plasmaResults/data
FASTQDUMP=/opt/apps/sratoolkit/2.3.4/bin/fastq-dump

for SRA in $FASTQPATH/*sra
do
    echo $FASTQDUMP --gzip --split-3 $SRA --outdir $FASTQPATH
done > command.sh

echo "#!/bin/bash
#
#$ -N extractSRA
#$ -pe 12way 12 
#$ -q normal
#$ -o extractSRA.o\$JOB_ID
#$ -l h_rt=96:00:00
#$ -V
#$ -cwd
#$ -A Exosome-RNA-seq
#------------------------------------------------------
module load launcher  gcc/4.7.1 bedtools sratoolkit
export EXECUTABLE=\$TACC_LAUNCHER_DIR/init_launcher 
export CONTROL_FILE=`pwd`/command.sh
export WORKDIR=$FASTQPATH

cd \$WORKDIR

\$TACC_LAUNCHER_DIR/paramrun \$EXECUTABLE \$CONTROL_FILE
" > job.sge
