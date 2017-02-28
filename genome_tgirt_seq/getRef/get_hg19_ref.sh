LINK=ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/chromFa.tar.gz
REF_DIR=/stor/work/Lambowitz/ref/hg19/Sequence/all_seq
OUT_TAR=$REF_DIR/chromFa.tar.gz
#wget --timestamping $LINK -O $OUT_TAR
#tar zxvf $OUT_TAR --directory $REF_DIR
cat $REF_DIR/*fa $REF_DIR/NC_007605.1.fasta > $REF_DIR/reference.fasta
rm $REF_DIR/*fa
bwa index $REF_DIR/reference.fasta
