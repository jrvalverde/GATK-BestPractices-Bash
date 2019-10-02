#!/bin/bash
#
# build indexes and interval lists
#
#	THIS SCRIPT MUST BE ADAPTED FOR EACH PROJECT CONSIDERING THE
#	PROJECT FILE NAMES.
#
#	IT MUST BE ADAPTED THOROUGHLY!!!! 
#
#	IT IS PROVIDED ONLY AS A REFERENCE EXAMPLE
#
#	Run this inside the project directory
#
#	It expects all sequence files to be in FASTQ format in a subdirectory
#	named 'fastq' and having long file names. It will generate a new
#	subdirectory 'fq' with shorter (yet still meaningful) file names
#	as symbolic links to the full-length-name files.
#

GATK3=~/contrib/gatk3
GATK4=~/contrib/gatk4
GATK_VERS=4

PRJ="ref"
# 
# cd $prj
# mkdir fq
# cd fq
# 
# for i in ../fastq/*.gz ; do 
# #    ln -s $i `basename $i \
# #              | sed -e "s/.*_S/${PRJ}_S/g" \
# #                    -e 's/_L00._/_/g' ` 
# 
# #    ln -s $i `basename $i | sed -e s/.*_17/ref/g -e s/\.fastq/_001\.fastq/g`
# done
# cd ..
#
 
#exit

cd data

if [ $GATK_VERS -eq 4 ] ; then
    if [ ! -e ucsc.hg19.dict ] ; then
        $GATK4/gatk CreateSequenceDictionary \
    	    -REFERENCE ucsc.hg19.fasta -OUTPUT ucsc.hg19.dict
    fi
#    if [ ! -e wg_hg19.dict ] ; then
#        $GATK4/gatk CreateSequenceDictionary \
#    	    -REFERENCE wg_hg19.fa -OUTPUT wg_hg19.dict
#    fi
else
    if [ ! -e ucsc.hg19.dict ] ; then
        java -Xms6000m -jar $GATK3/picard.jar CreateSequenceDictionary \
    	    REFERENCE=ucsc.hg19.fasta OUTPUT=ucsc.hg19.dict
    fi
#    if [ ! -e wg_hg19.dict ] ; then
#        java -Xms6000m -jar $GATK3/picard.jar CreateSequenceDictionary \
#    	    REFERENCE=wg_hg19.fa OUTPUT=wg_hg19.dict
#    fi
fi

if [ ! -e ucsc.hg19.fasta.fai ] ; then
    samtools faidx ucsc.hg19.fasta 
fi
if [ ! -e wg_hg19.fa.fai ] ; then
    samtools faidx wg_hg19.fa 
fi

if [ ! -e ucsc.hg19bwaidx.bwt ] ; then
    bwa index -a bwtsw -p ucsc.hg19bwaidx ucsc.hg19.fasta
fi

if [ ! -e ucsc.hg19.bt2.idx ] ; then
    bowtie-build ucsc.hg19.fasta ucsc.hg19.bt2.idx
fi

# lastdb ucsc.hg19.lst.idx ucsc.hg19.fasta

if [ ! -e wg_hg19bwaidx.bwt ] ; then
    bwa index -a bwtsw -p wg_hg19bwaidx wg_hg19.fasta
fi

if [ ! -e wg_hg19.bt2.idx ] ; then
    bowtie-build wg_hg19.fa wg_hg19.bt2.idx
fi

# lastdb wg_hg19.lst.idx wg_hg19.fa

cd -


# Prepare interval list
cd panel

if [ $GATK_VERS -eq 4 ] ; then
    $GATK4/gatk BedToIntervalList \
	    -I panel_65genes.bed \
            -O panel_65genes.ucsc.interval_list \
            -SD ../data/ucsc.hg19.fasta \
	    -R ../data/ucsc.hg19.fasta

    $GATK4/gatk BedToIntervalList \
	    -I panel_65genes.bed \
            -O panel_65genes.wg.interval_list \
            -SD ../data/wg_hg19.fa \
	    -R ../data/wg_hg19.fa
else
    java -Xms6000m -jar $GATK3/picard.jar BedToIntervalList \
	    INPUT=panel_65genes.bed \
	    OUTPUT=panel_65genes.ucsc.interval.list \
	    SEQUENCE_DICTIONARY=../data/ucsc.hg19.fasta

    java -Xms6000m -jar $GATK3/picard.jar BedToIntervalList \
	    INPUT=panel_65genes.bed \
	    OUTPUT=panel_65genes.wg.interval.list \
	    SEQUENCE_DICTIONARY=../data/wg_hg19.fa
fi

cd -

