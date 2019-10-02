#!/bin/bash
# gatk4-BP2019-PE-single-sample-GVFCs.sh

DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

## DATA DIRECTORY
data_dir="$HOME/work/lorena/data/"
gatk_bundle_dir="$data_dir/gatk-hg19-bundle"


## "PATHS",
bwa_path='/usr/bin/'
gatk_exec="$HOME/contrib/gatk4/gatk"
gotc_path="$HOME/contrib/gatk4/"
gitc_path="$HOME/contrib/gatk4/"
samtools='Needed if and when we add cramtobam()'
align_dir='align'
gvcf_dir='gvcf'

make_gvcf='YES'
contamination=0.0

java_opt="-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
command_mem_gb=16


# we need
#input_bam=$1
ref_name='ucsc.hg19'
#ref_name=$2
ref_dict="$data_dir/ucsc.hg19.dict"
ref_fasta="$data_dir/ucsc.hg19.fasta"

intervals_list=exome.intervals.list
#intervals_list=$3
wgs_calling_interval_list="exome.ucsc.hg19.interval_list"
wgs_evaluation_interval_list="exome.ucsc.hg19.interval_list"

dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"


# Identify GVFCs
function HaplotypeCallerGvcf_GATK4 {
    input_bam=$1
    ref_name=$2
    intervals_list=$3

    ref_fasta="$data_dir/$ref_name.fasta"
    out_dir='gvcf_bwa'
    mkdir -p $out_dir

    sample_basename=`basename $input_bam .aligned.duplicates_marked.recalibrated.bam`
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
#        -L ${interval_list} \
#        --dbsnp $gatk_bundle_dir/dbsnp_138.hg19.vcf \
    fi
}


# Validate a GVCF with -gvcf specific validation
function ValidateGVCF() {
    input_vcf=$1
  
    $CALL $gatk_exec --java-options "-Xms3000m" \
      ValidateVariants \
      -V ${input_vcf} \
      -R ${ref_fasta} \
      -L ${wgs_calling_interval_list} \
      -gvcf \
      --validation-type-to-exclude ALLELES \
      --dbsnp ${dbSNP_vcf}
}

# Collect variant calling metrics from GVCF output
function CollectGvcfCallingMetrics() {
    input_vcf=$1
    
    metrics_basename=`basename $1 .g.vcf.gz`
    output_basename=$gvcf_dir/$metrics_basename

    if [ ! -s ${output_basename}.variant_calling_summary_metrics \
      -o ! -s ${output_basenane}.variant_calling_detail_metrics ] ; then
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



function single_sample_variants() {
    input_vcf=$1

    # Validate the GVCF output of HaplotypeCaller
#    $CALL ValidateGVCF $input_vcf

    # QC the GVCF
    $CALL CollectGvcfCallingMetrics $input_vcf

}


for i in $align_dir/*.aligned.duplicates_marked.recalibrated.bam ; do
    $CALL HaplotypeCallerGvcf_GATK4 $i $ref_name $intervals_list
done

for i in $gvcf_dir/*.g.vcf.gz ; do
    $CALL single_sample_variants $i
done
