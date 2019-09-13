set -ue

if [ ! -e ucsc.hg19bwaidx.bwt ] ; then
    echo "Building UCSC HG 19 BWA index"
    bwa index -p ucsc.hg19bwaidx -a bwtsw ucsc.hg19.fasta
fi
if [ ! -e wg.hg19bwaidx.bwt ] ; then
    echo "Building WG HG 19 BWA index"
    bwa index -p wg.hg19bwaidx -a bwtsw wg_hg19.fa
fi


if [ ! -d bwa-align ] ; then mkdir bwa-align ; fi

for FILE in fastq-short/*_R1_* ; do
    echo Aligning with BWA: $FILE
    R1=$FILE
    R2=`echo $FILE | sed -e 's/_R1_/_R2_/g'`
    name=`basename $FILE`	# to remove the "fastq-short/" part
    OUT=`echo $name | sed -e 's/_R1_.*//g'`

    if [ ! -e bwa-align/${OUT}_R1_R2_ucsc-bwa.sam ] ; then
	REF=ucsc.hg19bwaidx
        bwa mem -M -t 6 $REF $R1 $R2 > bwa-align/${OUT}_R1_R2_ucsc-bwa.sam
    fi
    if [ ! -e bwa-align/${OUT}_R1_R2_wg-bwa.sam ] ; then
	REF=wg.hg19bwaidx
        bwa mem -M -t 6 $REF $R1 $R2 > bwa-align/${OUT}_R1_R2_wg-bwa.sam
    fi
done


