#!/bin/bash

# You should set the following four values appropriately
# they will be used to seed all the labels. They should 
# also be consistent with the configuration file used later
# to apply the GATK workflow.

prefix="PR1"
library="PRORADIUM1"
sequencer="INH"
date='2018-01-01T00:00:00+0100'

#prefix="PR2"
#lirary='PRORADIUM2'
#sequencer='OUT'		# outside
#date='2019-01-01T00:00:00+0100'



GATK3=~/contrib/gatk3
GATK4=~/contrib/gatk4

R1=$1
R2=`echo $1 | sed -e 's/_R1_/_R2_/g'`
F=`basename $R1 .gz`	# remove .gz from name just in case it is compressed
BAM=../uBam/`basename $F .fastq | sed -e 's/_R1_/_R1+2_/g'`.bam

sample_id=`echo $F | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
#READ_GROUP_NAME=$prefix-${sample_id}
READ_GROUP_NAME=$prefix${sample_id}
#SAMPLE_NAME=${sample_id}-`echo $F | cut -d'-' -f1`
SAMPLE_NAME=${sample_id}`echo $F | cut -d'-' -f1`
platform_unit=`zcat $R1 | head -1 | cut -d':' -f1 | tr -d '@'`
### Looks like GATK is using PLATFORM_UNIT for the read group name!!!
### at least on the data-processing workflow BQSR reports!!!
#PLATFORM_UNIT=$platform_unit
PLATFORM_UNIT=${sample_id}${platform_unit}
#
LIBRARY_NAME=$library
SEQUENCING_CENTER=$sequencer			# in-house
RUN_DATE=$date

mkdir -p ../uBam/

echo java -Xmx8G -jar $GATK3/picard.jar FastqToSam \
    FASTQ=$R1 \
    FASTQ2=$R2 \
    OUTPUT=$BAM \
    READ_GROUP_NAME=$READ_GROUP_NAME \
    SAMPLE_NAME=$SAMPLE_NAME \
    LIBRARY_NAME=$LIBRARY_NAME \
    PLATFORM=illumina \
    SEQUENCING_CENTER=$SEQUENCING_CENTER \
    RUN_DATE=$RUN_DATE \
    PLATFORM_UNIT=$PLATFORM_UNIT 

$GATK4/gatk FastqToSam \
    --FASTQ  $R1 \
    --FASTQ2 $R2 \
    --OUTPUT $BAM \
    --READ_GROUP_NAME   $READ_GROUP_NAME \
    --SAMPLE_NAME       $SAMPLE_NAME \
    --LIBRARY_NAME      $LIBRARY_NAME \
    --PLATFORM          illumina \
    --SEQUENCING_CENTER $SEQUENCING_CENTER \
    --RUN_DATE          "$RUN_DATE" \
    --PLATFORM_UNIT     $PLATFORM_UNIT 
