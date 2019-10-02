## Copyright Broad Institute, 2019
## 
## The haplotypecaller-gvcf-gatk4 workflow runs the HaplotypeCaller tool
## from GATK4 in GVCF mode on a single sample according to GATK Best Practices.
## When executed the workflow scatters the HaplotypeCaller tool over a sample
## using an intervals list file. The output file produced will be a
## single gvcf file which can be used by the joint-discovery workflow.
##
## Requirements/expectations :
## - One analysis-ready BAM file for a single sample (as identified in RG:SM)
## - Set of variant calling intervals lists for the scatter, provided in a file
##
## Outputs :
## - One GVCF file and its index
##
## Cromwell version support 
## - Successfully tested on v37
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation.
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

# we need
#input_bam=$1
ref_name='ucsc.hg19'
#ref_name=$2
intervals_list=exome.intervals.list
#intervals_list=$3


## DATA DIRECTORY
data_dir="$HOME/work/lorena/data/"
gatk_bundle_dir="$data_dir/gatk-hg19-bundle"


## "PATHS",
bwa_path='/usr/bin/'
gatk_exec="$HOME/contrib/gatk4/gatk"
gotc_path="$HOME/contrib/gatk4/"
gitc_path="$HOME/contrib/gatk4/"
samtools='Needed if and when we add cramtobam()'

make_gvcf='YES'
contamination=0.0

java_opt="-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
command_mem_gb=16

# WORKFLOW DEFINITION 
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

for i in align_bwa/*.aligned.duplicates_marked.recalibrated.bam ; do
    $CALL HaplotypeCallerGvcf_GATK4 $i $ref_name $intervals_list
done
