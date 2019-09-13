#!/bin/bash
#
#	Check quality using fastQC
#
# (C) Lorena Magraner, 2019
# License: GPL
#
# This script assumes that there is a subfolder named 'fastq' containing
# all the compressed fastq files.
#
# We will generate a subfolder named 'qc' containing the quality ZIP
# files. 
#
# Requires the AdapterRemoval software installed

# this is where AdapterRemoval is installed
alias AdapterRemoval=~/contrib/AdapterRemovalV2/bin/AdapterRemoval

set -ue

if [ ! -d qc ] ; then mkdir qc ; fi
if [ ! -d qc/adapter ] ; then mkdir qc/adapter ; fi

cd qc

for FILE in ../fq/*_R1_* ; do
    echo FastQC Quality Control for  $FILE
    R1=$FILE
    R1out=`basename $R1 .fastq.gz`_fastqc
    R2=`echo $FILE | sed -e 's/_R1_/_R2_/g'`
    R2out=`basename $R2 .fastq.gz`_fastqc

    if [ ! -d $R1out ] ; then
        fastqc --noextract "$R1" --outdir .
        unzip $R1out.zip
        rm $R1out.zip
    fi

    if [ ! -d $R2out ] ; then
        fastqc --noextract "$R2" --outdir .
        unzip $R2out.zip
        rm $R2out.zip
    fi

    # Identify potential adapters
    ADAPT=adapter/`basename $R1 .fastq.gz`_adapters.txt
    if [ ! -d adapter ] ; then mkdir adapter ; fi
    if [ ! -e $ADAPT ] ; then
        ~/contrib/AdapterRemovalV2/bin/AdapterRemoval \
            --identify-adapters --file1 $R1 --file2 $R2 \
            |& tee $ADAPT 2>&1
    fi
    
done

# Make a global report
if [ ! -e multiqc_report.html ] ; then
    multiqc .
fi

# Check adapters
cd adapter

echo "adapter-D" > ../adapters.txt
grep "adapter1:" * | cut -d' ' -f5 | uniq >> ../adapters.txt
echo "adapter-R" >> ../adapters.txt
grep "adapter2:" * | cut -d' ' -f5 | uniq >> ../adapters.txt

echo "9bp-5'kmer-D-1st" >> ../adapters.txt
grep " 1:" * | cut -d':' -f3 | cut -d'=' -f1 | sort | uniq >> ../adapters.txt
echo "9bp-5'kmer-R-1st" >> ../adapters.txt
grep " 2:" * | cut -d':' -f3 | cut -d'=' -f1 | sort | uniq >> ../adapters.txt
