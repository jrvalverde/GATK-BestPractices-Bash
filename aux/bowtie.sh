#!/bin/bash
ref=$1
s1=$2
s2=$3

set -e

db=`basename $ref .fas`

if [ ! -e $ref.fai ] ; then
    bowtie-build $ref $db
fi

if [ ! -e EPEC.bowtie.map ] ; then
    bowtie -t $db -1 $s1 -2 $s2 EPEC.bowtie.map
fi
if [ ! -e EPEC.bowtie.sam ] ; then
    bowtie -t -p 8 -S -q -n 3 -e 90 -l 50 -y $db -1 $s1 -2 $s2 EPEC.bowtie.sam
fi

#samtools view -bS -q 99 EPEC.bowtie.sam > EPEC.bowtie.bam
if [ ! -e EPEC.bowtie.bam ] ; then
    samtools import $ref EPEC.bowtie.sam EPEC.bowtie.bam
fi
if [ ! -e EPEC.bowtie.sorted.bam ] ; then
    samtools sort EPEC.bowtie.bam EPEC.bowtie.sorted
fi
if [ ! -e EPEC.bowtie.sorted.bam.bai ] ; then
    samtools index EPEC.bowtie.sorted.bam
fi
if [ ! -e EPEC.bowtie.pileup ] ; then
    samtools pileup -cf $ref EPEC.bowtie.sorted.bam > EPEC.bowtie.pileup
fi
if [ ! -e EPEC.bowtie.pileup.fastq ] ; then
    perl /usr/share/samtools/samtools.pl pileup2fq EPEC.bowtie.pileup > EPEC.bowtie.pileup.fastq
fi

if [ ! -e EPEC.bowtie.sorted.bcf ] ; then
    #samtools mpileup -u -f $ref EPEC.bowtie.sorted.bam > EPEC.bowtie.sorted.bcf
    samtools pileup -f $ref EPEC.bowtie.sorted.bam > EPEC.bowtie.sorted.bcf
fi

if [ ! -e EPEC.bowtie.sorted.vcf ] ; then
    bcftools view -v -c -g EPEC.bowtie.sorted.bcf > EPEC.bowtie.sorted.vcf
fi
if [ ! -e EPEC.bowtie.sorted.fastq ] ; then
    perl /usr/share/samtools/vcfutils.pl vcf2fq EPEC.bowtie.sorted.vcf > EPEC.bowtie.sorted.fastq
fi
