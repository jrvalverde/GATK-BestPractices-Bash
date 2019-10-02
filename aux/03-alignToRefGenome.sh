#!/bin/bash
#
#	Align sequences to the reference genome(s) using BWA
#
# This script is to be run inside the project directory. It expects
# experiment sequences to be stored as fastq files in a subdirectory
# named "fq" and will store the results in a subfolder named '$align_dir'
# as SAM files
#

set -ue

USE_BWA="YES"
USE_BOWTIE2="NO"

align_dir=align.gatk

if [ "$USE_BWA" == "YES" ] ; then

    if [ ! -d $align_dir ] ; then mkdir $align_dir ; fi

    for FILE in fq/*_R1_* ; do
        echo Aligning with BWA: $FILE
        R1=$FILE
        R2=`echo $FILE | sed -e 's/_R1_/_R2_/g'`
        name=`basename $FILE`	# to remove the "fastq-short/" part
        OUT=`echo $name | sed -e 's/_R1_.*//g'`

	sample=`basename $FILE | sed -e 's/^*_S//g' -e 's/_.*//g'`
	# This should be used for alignment
        rginfo="@RG\tID:1\tSM:$sample\tPL:illumina\tLB:proradium\tPU:unit1"

        if [ ! -e $align_dir/${OUT}_R1_R2-ucsc-bwa.sam ] ; then
	    REF=../data/ucsc.hg19bwaidx
	    bwa mem \
                -R "$rginfo" \
                -M -t 6 $REF $R1 $R2 \
                > $align_dir/${OUT}_R1_R2-ucsc-bwa.sam 
        fi

#        if [ ! -e $align_dir/${OUT}_R1_R2-wg-bwa.sam ] ; then
	if [ "YES" == "NO" ] ; then
	    REF=../data/wg.hg19bwaidx
            bwa mem \
            	-R "$rginfo" \
            	-M -t 6 $REF $R1 $R2 \
                > $align_dir/${OUT}_R1_R2-wg-bwa.sam
        fi


    done
fi


if [ "$USE_BOWTIE2" == "YES" ] ; then
    # ALIGN WITH BOWTIE2

    if [ ! -d align-bwt2 ] ; then mkdir align-bwt2 ; fi

    for FILE in fq/*_R1_* ; do
        echo Aligning with BOWTIE2: $FILE
        R1=$FILE
        R2=`echo $FILE | sed -e 's/_R1_/_R2_/g'`
        name=`basename $FILE`	# to remove the "fastq-short/" part
        OUT=`echo $name | sed -e 's/_R1_.*//g'`

        if [ ! -e align-bwt2/${OUT}_R1_R2-ucsc-bwt2.sam ] ; then
	    REF=../data/ucsc.hg19.bt2.idx
            bowtie2 -p 8 --qc-filter -x $REF \
	        -1 $R1 -2 $R2 \
                -S align-bwt2/${OUT}_R1_R2-ucsc-bwt2.sam
        fi
        
        if [ "YES" == "NO" ] ; then
        #if [ ! -e align-bwt2/${OUT}_R1_R2-wg-bwt2.sam ] ; then
	    REF=../data/wg_hg19.bt2.idx
            bowtie2 -p 8 --qc-filter -x $REF \
	        -1 $R1 -2 $R2 \
                -S align-bwt2/${OUT}_R1_R2-wg-bwt2.sam
        fi
    done
fi
