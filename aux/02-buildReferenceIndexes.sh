#!/bin/bash
#
#	Create BWA indexes for the reference genome sequence
#
#	This script should be run inside the directory containing
# the reference genome


if [ ! -e ucsc.hg19bwaidx.bwt ] ; then
    echo "Building UCSC HG 19 BWA index"
    bwa index -p ucsc.hg19bwaidx -a bwtsw ucsc.hg19.fasta
fi
