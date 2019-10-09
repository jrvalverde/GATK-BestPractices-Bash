#!/bin/bash
#
# This script has been developed at Scientific Computing Service of
# the National Biotechnology Center of the Spanish Research Council 
# (CNB-CSIC) and it is therefore the property of CNB-CSIC.
#
# It is based on Cromwell scripts developed by the Broad Institute.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
#---------------------------------------------------------------------
## Copyright CNB-CSIC, 2019
## Copyright Broad Institute, 2019
#

# NOTE: the actual workflow is the last function defined in this file.
#	THIS WORKFLOW IS NOW SUPERSEDED BY 
#		gatk4-BP2019-pe-single-sample-wf+filtering.sh

set -x

DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

################### THESE PROVIDE DEFAULT VALUES ###################
################## OVERRIDEN IN gatk4-BP2019.in.sh #################
## DATA DIRECTORY
data_dir="./hg19"
gatk_bundle_dir="$data_dir/gatk-bundle"
## "PATHS",
bwa_path='/usr/bin/'
gatk_exec="$HOME/contrib/gatk4/gatk"
gotc_path="$HOME/contrib/gatk4/"
gitc_path="$HOME/contrib/gatk3/"
samtools='Needed if and when we add cramtobam()'
align_dir='align.gatk'
gvcf_dir='g.vcf'
# "AUXILIARY FILES"
ref_name='ucsc.hg19'
study_interval_list=".$panel_dir/panel5-120.ucsc.hg19.interval_list"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"

# get actual preferences from a local file
if [ -s gatk4-BP2019.in.sh ] ; then
    . gatk4-BP2019.in.sh
else
    cp `dirname $0`/gatk4-BP2019.in.sh \
       ./gatk4-BP2019.in.sh
    echo "I have created a preferences file in this directory:"
    echo ""
    echo "        gatk4-BP2019.in.sh"
    echo ""
    echo "Please, open it in a text editor, adapt it to your needs and"
    echo "run this command again."
    exit
fi

make_gvcf='YES'
contamination=0.0

    if [ "$make_gvcf" == 'YES' ] ; then
        gvcf_dir="g.vcf"
    else
        gvcf_dir="vcf"
    fi

java_opt="-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
command_mem_gb=16

# we need
ref_dict="$data_dir/$ref_name.dict"
ref_fasta="$data_dir/$ref_name.fasta"

wgs_calling_interval_list="$study_interval_list"
wgs_evaluation_interval_list="$study_interval_list"


## HaplotypeCallerGvcf_GATK4()
#
# Usage:
# 	HaplotypeCallerGvcf_GATK4 recal.bam ref.name interval_list
#
# Identify GVFCs (at one point we also supported GATK V3, and even V2, 
# but we no longer do)
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function HaplotypeCallerGvcf_GATK4() {
    if [ $# -ne 3 ] ; then
        echo "Usage: ${FUNCNAME[0]} recal.bam ref.name interval_list"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3"
    fi
    local input_bam=$1
    local ref_name=$2
    local interval_list=$3

    local ref_fasta="$data_dir/$ref_name.fasta"
    local out_dir='g.vcf'
    mkdir -p $out_dir

    sample_basename=`basename $input_bam .bam`
    vcf_basename=$sample_basename
    if [ "$make_gvcf" == 'YES' ] ; then
        output_suffix="g.vcf.gz"
        gvcf_opt="-ERC GVCF"
    else
        output_suffix="vcf.gz"
        gvcf_opt=""
    fi
    output_filename=${out_dir}/${vcf_basename}.${output_suffix}

    # Generate GVCF
    set -e

    if [ ! -s ${output_filename} ] ; then
    $CALL ${gatk_exec} --java-options "-DGATK_STACKTRACE_ON_USER_EXCEPTION=true -Xmx${command_mem_gb}G ${java_opt}" \
        HaplotypeCaller \
        -R ${ref_fasta} \
        -I ${input_bam} \
        -O ${output_filename} \
	$gvcf_opt
#
# Additional possible options (the use has been changing quite a lot and there
# are tons of variations reported as standard workflow in the Internet). We are
# using just the ones listed in the latest (as of 2019) workflows.
#        -L ${interval_list} \
#        --dbsnp $gatk_bundle_dir/dbsnp_138.hg19.vcf :
#		rsIDs from this file are used to populate the ID column of 
#		the output. Also, the DB INFO flag will be set when a
#		ppropriate. dbSNP is not used in any way for the calculations 
#		themselves.
#	--genotyping_mode DISCOVERY 
#	--doNotRunPhysicalPhasing 
#	--variant_index_type LINEAR 
#	--variant_index_parameter 128000



    fi
}


## ValidateGVCF()
#
# Usage:
#	ValidateGVCF input_vcf
#
# Validate a GVCF with gvcf specific validation
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function ValidateGVCF() {
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} input_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1"
    fi
    local input_vcf=$1
    
    echo ">>> Validating $input_vcf"
    #return   #<<<########################### uncomment to avoid this step
    # NOTE: this may fail validation on exome data
    # which is why we use --warnOnErrors
    $CALL $gatk_exec --java-options "-Xms3000m" \
        ValidateVariants \
        -V ${input_vcf} \
        -R ${ref_fasta} \
        -L ${wgs_calling_interval_list} \
        -gvcf \
        --validation-type-to-exclude ALLELES \
        --dbsnp ${dbSNP_vcf} \
        --warnOnErrors
}


## CollectGvcfCallingMetrics()
#
# Usage:
#	CollectGvcfCallingMetrics input.vcf
#
# Collect variant calling metrics from GVCF output
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function CollectGvcfCallingMetrics() {
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} input_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1"
    fi
    local input_vcf=$1
    
    local metrics_basename=`basename $1 .g.vcf.gz`
    local output_basename=$gvcf_dir/$metrics_basename

    if [ ! -s ${output_basename}.variant_calling_summary_metrics \
      -o ! -s ${output_basename}.variant_calling_detail_metrics ] ; then
        # For GATK v3
        #$CALL java -Xms2000m -jar ${gitc_path}/picard.jar \
        #  CollectVariantCallingMetrics \
        #  INPUT=${input_vcf} \
        #  OUTPUT=${output_basename} \
        #  DBSNP=${dbSNP_vcf} \
        #  SEQUENCE_DICTIONARY=${ref_dict} \
        #  TARGET_INTERVALS=${wgs_evaluation_interval_list} \
        #  GVCF_INPUT=true

	# for GATK v4
        $CALL $gatk_exec --java-options "-Xms2000m" \
            CollectVariantCallingMetrics \
            --INPUT $input_vcf \
            --OUTPUT $output_basename \
            --DBSNP $dbSNP_vcf \
            --TARGET_INTERVALS ${wgs_evaluation_interval_list} \
            --SEQUENCE_DICTIONARY $ref_dict \
            --GVCF_INPUT true
    fi

}


## AnnotateVariants()
#
# Usage:
#	AnnotateVariants input_bam input_vcf
#
#	Annotate the variants with the help of the input BAM file and
# dbSNP database. Note that we could also run snpEff and use the vcf file
# produced to include additional annotation in this function. We should
# actually consider it at some point.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function AnnotateVariants() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} input_bam input_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    # in principle this should not be needed
    local input_bam=$1
    local input_vcf=$2
    # substitute if match is at end of string
    local output=${input_vcf/%$output_suffix/ann.g.vcf}
    
    if [ ! -s $output ] ; then
        # Note that the majority of these annotations will not be used!!!
        # `# comment` is just a way to insert comments in line.
        $CALL $gatk_exec \
    	    VariantAnnotator \
            -R $ref_fasta \
            -I $input_bam \
            -V $input_vcf \
            -O $output \
            -A AlleleFraction \
            -A BaseQuality \
            -A BaseQualityRankSumTest \
            -A ChromosomeCounts \
            -A CountNs \
            -A Coverage \
            -A DepthPerAlleleBySample \
            -A DepthPerSampleHC \
            -A FisherStrand \
            -A FragmentLength \
            -A GenotypeSummaries \
            -A QualByDepth \
            -A ReadPosRankSumTest \
            -A MappingQuality \
            -A MappingQualityRankSumTest \
            -A MappingQualityZero \
            -A AS_BaseQualityRankSumTest \
            -A AS_FisherStrand \
            -A AS_InbreedingCoeff \
            -A AS_MappingQualityRankSumTest \
            -A AS_QualByDepth \
            `# -A AS_RMSMappingQuality` \
            -A AS_ReadPosRankSumTest \
            -A AS_StrandOddsRatio \
            -L $input_vcf \
            --dbsnp $dbSNP_vcf
    fi
}

## single_sample_gvcfs()
#
# Usage:
#	single_sample_gvcfs recal.bam ref.name interval_list
#
#	This is the actual workflow. Produce GVCFs from the preprocessed,
# recalibrated BAM aignment files.
#
#	THIS WORKFLOW IS NOW SUPERSEDED BY 
#		gatk4-BP2019-pe-single-sample-wf+filtering.sh
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function single_sample_gvcfs() {
    if [ $# -ne 3 ] ; then
        echo "Usage: ${FUNCNAME[0]} recal.bam ref.name interval_list"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3"
    fi
    local ali_srt_dedup_recal_bam=$1
    local ref_name=$2
    local intervals_list=$3
    local vcf=''
    if [ "$make_gvcf" == 'YES' ] ; then
        output_suffix='g.vcf.gz'
    else
        output_suffix='vcf.gz'
    fi
    vcf=$gvcf_dir/`basename $ali_srt_dedup_recal_bam bam`$output_suffix
    

    $CALL HaplotypeCallerGvcf_GATK4 $ali_srt_dedup_recal_bam \
    	$ref_name \
        $intervals_list
    echo ">>> $vcf"
    $CALL ValidateGVCF "$vfc"
    $CALL CollectGvcfCallingMetrics "$vcf"
    $CALL AnnotateVariants "$ali_srt_dedup_recal_bam" "$vcf"

}


#########################################################################
#### Do the actual work

# if we are called with no arguments process all recalibrated files
if [ $# -eq 0 ] ; then
    #$0 $align_dir/*.recal*.bam
    # we cannot pass multi-sample files to HaplotypeCaller
    $0 `ls $align_dir/*.recal*.bam | grep -v "PR"`
    exit
elif [ $# -gt 1 ] ; then
    # if we get more than one BAM file, process each one by one
    for i in  "$@" ; do
        $0 $i
    done
    exit
fi
# $# must be exactly one BAM file to process, do it
recalibrated_bam=$1

echo ">>> $0 processing $1"

single_sample_gvcfs $recalibrated_bam $ref_name $study_interval_list

#########################################################################
