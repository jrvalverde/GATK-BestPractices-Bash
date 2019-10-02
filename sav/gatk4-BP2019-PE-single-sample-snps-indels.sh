#!/bin/bash


DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

data_dir="$HOME/work/lorena/data/"
gatk_bundle_dir="$data_dir/gatk-hg19-bundle"
gatk_exec="$HOME/contrib/gatk4/gatk"
gotc_path="$HOME/contrib/gatk4/"
gitc_path="$HOME/contrib/gatk4/"
gvcf_dir='gvcf_bwa'

ref_dict="$data_dir/ucsc.hg19.dict"
ref_fasta="$data_dir/ucsc.hg19.fasta"
wgs_calling_interval_list="exome.ucsc.hg19.interval_list"
wgs_evaluation_interval_list="exome.ucsc.hg19.interval_list"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"


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
        #$CALL java -Xms2000m -jar ${gitc_path}/picard.jar \
        #  CollectVariantCallingMetrics \
        #  INPUT=${input_vcf} \
        #  OUTPUT=${output_basename} \
        #  DBSNP=${dbSNP_vcf} \
        #  SEQUENCE_DICTIONARY=${ref_dict} \
        #  TARGET_INTERVALS=${wgs_evaluation_interval_list} \
        #  GVCF_INPUT=true

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

for i in $gvcf_dir/*.g.vcf.gz ; do
    single_sample_variants $i
done
