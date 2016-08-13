#/usr/bin/bash
# copy from https://github.com/chapmanb/bcbio-nextgen/blob/master/config/examples/NA12878-trio-wgs-validate-getdata.sh


DATAPATH=/stor/work/Lambowitz/Data/NGS/NA12878_gDNA

curl -o $DATAPATH/NA12878_1.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_1.fastq.gz
curl -o $DATAPATH/NA12878_2.fastq.gz ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR194/ERR194147/ERR194147_2.fastq.gz

