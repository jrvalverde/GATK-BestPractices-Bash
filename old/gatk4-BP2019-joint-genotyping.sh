#!/bin/bash

data_dir="/home/scientific/work/lorena/data/"
gatk_bundle_dir="$data_dir/gatk-hg19-bundle"
gatk_exec="/home/scientific/contrib/gatk4/gatk"
gatk_path="/home/scientific/contrib/gatk4/"
gotc_path="/home/scientific/contrib/gatk4/"
gitc_path="/home/scientific/contrib/gatk4/"
align_dir='align'
gvcf_dir='gvcf'
workspace_dir_name='joint_gvcf'
out_dir='joint_gvcf_analysis'

ref_name=ucsc.hg19
ref_fasta="$data_dir/ucsc.hg19.fasta"
ref_dict="$data_dir/ucsc.hg19.dict"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"
hapmap_vcf="$gatk_bundle_dir/hapmap_3.3.hg19.sites.vcf"
omni_vcf="$gatk_bundle_dir/1000G_omni2.5.hg19.sites.vcf"
Mills_vcf="$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
G1000_vcf="$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"

SNP_VQSR_downsampleFactor=10
snp_filter_level=99.7
indel_filter_level=99.7
indel_recalibration_annotation_values=("FS" "ReadPosRankSum" "MQRankSum" "QD" "SOR" "DP")
snp_recalibration_annotation_values=("QD" "MQRankSum" "ReadPosRankSum" "FS" "MQ" "SOR" "DP")

uBam_list='uBam.list'
picard_interval_list='exome.ucsc.hg19.interval_list'


# prepare sample name map from uBam.list
if [ ! -e sample_name.map ] ; then
    cat $uBam_list | while read bam ; do
        F=`basename $bam .bam`	# remove bam extension
        sample_id=`echo $F | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
        #SAMPLE_NAME=${sample_id}-`echo $F | cut -d'-' -f1`
        SAMPLE_NAME=${sample_id}`echo $F | cut -d'-' -f1`
        file_name=$gvcf_dir/$F.g.vcf.gz
        echo "$SAMPLE_NAME	$file_name"
    done > sample_name.map
fi
sample_name_map='sample_name.map'

# prepare arguments with gvcf list
gvcf_arg=''
cat $sample_name_map | while read sample file ; do
    gvcf_arg="${gvcf_arg} -V $file "
done


# Step 1: import gvcfs
intervals=$picard_interval_list

if [ ! -s ${workspace_dir_name}.tar ] ; then
    echo "Importing GVCFs"
    set -e
    
    rm -rf ${workspace_dir_name}

    batch_size=0	# read all samples at once
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    $gatk_exec --java-options "-Xmx24g -Xms24g" \
    GenomicsDBImport \
    --genomicsdb-workspace-path ${workspace_dir_name} \
    --batch-size ${batch_size} \
    -L ${intervals} \
    --sample-name-map ${sample_name_map} \
    --reader-threads 5 \
    --interval-padding 500 \
    --reference $ref_fasta \
    --validate-sample-name-map

    tar -cf ${workspace_dir_name}.tar ${workspace_dir_name}

fi


# Step 2: Genotype GVCFs
WORKSPACE=${workspace_dir_name}
joint_vcf="$out_dir/joint_gvcf.vcf.gz"

if [ ! -s $joint_vcf ] ; then
    echo "Genotyping GVCFs"
    set -e

    $gatk_exec --java-options "-Xmx5g -Xms5g" \
     GenotypeGVCFs \
         -R ${ref_fasta} \
         -O ${joint_vcf} \
         -D ${dbSNP_vcf} \
         --only-output-calls-starting-in-intervals \
         --use-new-qual-calculator \
         -V gendb://$WORKSPACE \
         -L ${intervals}
         #-G StandardAnnotation \ Unneeded (it is the default)
fi

if [ ! -s $out_dir/joint_gvcf_combined.g.vcf ] ; then
    # generate a flat multisample GVCF file from the GenomicsDB
    $gatk_exec SelectVariants \
        -R $ref_fasta \
        -V gendb://$WORKSPACE \
        -O $out_dir/joint_gvcf_combined.g.vcf

fi


excess_het_thresold=54.69
variant_filtered_vcf="$out_dir/joint_gvcf.variant_filtered.vcf.gz"
variant_filtered_nocall_vcf="$out_dir/joint_gvcf.variant_filtered.nocall.vcf.gz"
sites_only_vcf="$out_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz"
sites_only_nocall_vcf="$out_dir/joint_gvcf.variant_filtered.sites_only.nocall.vcf.gz"

# Hard-filter a large cohort callset on ExcessHet using VariantFiltration
if [ ! -s ${variant_filtered_vcf} ] ; then
    echo "Hard filtering"
    $gatk_exec --java-options "-Xmx3g -Xms3g" \
        VariantFiltration \
        --filter-expression "ExcessHet > 54.69" \
        --filter-name ExcessHet \
        -O ${variant_filtered_vcf} \
        -V ${joint_vcf}
#        --filter-expression "ExcessHet > ${excess_het_threshold}" \
fi

# transform filtered genotypes to no-call
if [ ! -s $variant_filtered_nocall_vcf ] ; then
    echo "Making nocall file"
    $gatk_exec SelectVariants \
    -V $variant_filtered_vcf \
    --set-filtered-gt-to-nocall \
    -O $variant_filtered_nocall_vcf
fi

# Create sites-only VCF with MakeSitesOnlyVcf
if [ ! -s ${sites_only_vcf} ] ; then
    echo "Making sites-only VCF"
    #java -Xmx3g -Xms3g -jar $gitc_path/picard.jar \
    $gatk_exec --java-options "-Xmx3g -Xms3g" \
        MakeSitesOnlyVcf \
        -INPUT ${variant_filtered_vcf} \
        -OUTPUT ${sites_only_vcf}
fi

if [ ! -s ${sites_only_nocall_vcf} ] ; then
    echo "Making no-call sites-only VCF"
    #java -Xmx3g -Xms3g -jar $gitc_path/picard.jar \
    $gatk_exec --java-options "-Xmx3g -Xms3g" \
        MakeSitesOnlyVcf \
        -INPUT ${variant_filtered_nocall_vcf} \
        -OUTPUT ${sites_only_nocall_vcf}
fi

# Calculate INDELs using VariantRecalibrator
indels_recal="$out_dir/joint_gvfcs.indels.recal.vcf"
indels_model_report="$out_dir/joint_gvfcs.indels.model.report"
indels_tranches="$out_dir/joint_gvfcs.indels.tranches"

if [ ! -s $indels_recal -o ! -s $indels_tranches ] ; then
    echo "Recalibrating INDELS"
    $gatk_exec --java-options "-Xmx24g -Xms24g" \
        VariantRecalibrator \
        -V ${sites_only_vcf} \
        -O ${indels_recal} \
        --tranches-file ${indels_tranches} \
        --trust-all-polymorphic \
        -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
        -an FS -an ReadPosRankSum -an QD -an SOR -an DP \
        -mode INDEL \
        --max-gaussians 4 \
        -resource:mills,known=false,training=true,truth=true,prior=12 \
	    ${Mills_vcf} \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2 \
	    ${dbSNP_vcf} \
	-L $intervals

#-an "FS" -ab "ReadPosRankSum" -an "MQRankSum" -an "QD" -an "SOR" -an "DP" \
#        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 \
#	    ${axiomPoly_vcf} \
fi
if [ ! -s $indels_recal -o ! -s $indels_tranches ] ; then
        echo "Recalibrating INDELS AGAIN!"
        $gatk_exec --java-options "-Xmx3g -Xms3g" \
            VariantRecalibrator \
            -V ${sites_only_vcf} \
            -O ${indels_recal} \
            --tranches-file ${indels_tranches} \
            --trust-all-polymorphic \
            -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
            -an FS -an ReadPosRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            -L $intervals \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15 \
        	${hapmap_vcf} \
            -resource:mills,known=false,training=true,truth=true,prior=12 \
                  ${Mills_vcf} \
            -resource:omni,known=false,training=true,truth=true,prior=12 \
        	${omni_vcf} \
            -resource:1000G,known=false,training=true,truth=false,prior=10 \
        	${G1000_vcf} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        	${dbSNP_vcf} 

fi

#  Calculate VQSLOD tranches for SNPs using VariantRecalibrator 
snps_model_report="$out_dir/joint_gvfcs.snps.model.report"
snps_recal="$out_dir/joint_gvcfs.snps.recal.vcf"
snps_tranches="$out_dir/joint_gvfcs.snps.tranches"
downsampleFactor=$SNP_VQSR_downsampleFactor

num_gvcfs=`wc -l $sample_name_map | cut -d' ' -f1`

if [ $num_gvcfs -gt 1000 ] ; then
    if [ ! -s $snps_model_report -o ! -s $snps_tranches ] ; then
    echo "Creating model for SNPs"
    $gatk_exec --java-options "-Xmx24g -Xms24g" \
        VariantRecalibrator \
        -V ${sites_only_vcf} \
        -O $snps_recal \
        --tranches-file $snps_tranches \
        --trust-all-polymorphic \
        -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
        -an FS -an ReadPosRankSum -an QD -an SOR -an DP \
        -mode SNP \
        --sample-every-Nth-variant ${downsampleFactor} \
        --output-model ${snps_model_report} \
        --max-gaussians 6 \
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
            $hapmap_vcf \
        -resource:mills,known=false,training=true,truth=true,prior=12.0 \
            ${Mills_vcf} \
        -resource:omni,known=false,training=true,truth=false,prior=12.0 \
            $omni_vcf \
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
            $G1000_vcf \
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
            $dbSNP_vcf \
	-L $intervals \
        #-an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
      
    fi
    recalibration_filename=joint_gvfc.snp.recal
    tranches_filename=joint_gvfc.snp.tranches
    
    indels_recal="joint_gvcfs_indels.recal"
    indels_tranches="joint_gvcfs_indels.tranches"

    # Recalibrate SNPs with  the previous model
    $gatk_exec --java-options "-Xmx3g -Xms3g" \
      VariantRecalibrator \
      -V ${sites_only_nocall_vcf} \
      -O ${snps_recal} \
      --tranches-file ${snps_tranches} \
      --trust-all-polymorphic \
      -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
      -an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
      -mode SNP \
      --input-model $snps_model_report \
      --max-gaussians 6 \
      -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
          $hapmap_vcf \
      -resource:mills,known=false,training=true,truth=true,prior=12.0 \
          ${Mills_vcf} \
      -resource:omni,known=false,training=true,truth=false,prior=12.0 \
          $omni_vcf \
      -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
          $G1000_vcf \
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
          $dbSNP_vcf 

else
    if [ ! -e "$snps_recal" -o ! -e "$snps_tranches" ] ; then
	echo "Recalibrating SNPs directly"
	# Recalibrating SNPs directly without model
        $gatk_exec --java-options "-Xmx3g -Xms3g" \
            VariantRecalibrator \
            -V ${sites_only_vcf} \
            -O ${snps_recal} \
            --tranches-file ${snps_tranches} \
            --trust-all-polymorphic \
            -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
            -an FS -an ReadPosRankSum -an QD -an SOR -an DP \
            -mode SNP \
            -L $intervals \
            --max-gaussians 6 \
            -resource:hapmap,known=false,training=true,truth=true,prior=15 \
        	${hapmap_vcf} \
            -resource:mills,known=false,training=true,truth=true,prior=12 \
                  ${Mills_vcf} \
            -resource:omni,known=false,training=true,truth=true,prior=12 \
        	${omni_vcf} \
            -resource:1000G,known=false,training=true,truth=false,prior=10 \
        	${G1000_vcf} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        	${dbSNP_vcf} 
    fi
#-an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
# MQ has no data and so MQRankSum has zero variance
fi


# Apply Recalibration

if [ -e indels_recal -a ! -s tmp.indel.recalibrated.vcf ] ; then
    set -e
    echo "Applying INDEL recalibration"

    $gatk_exec --java-options "-Xmx5g -Xms5g" \
	ApplyVQSR \
	-O $out_dir/tmp.indel.recalibrated.vcf \
	-V ${sites_only_vcf} \
	--recal-file ${indels_recal} \
	--tranches-file ${indels_tranches} \
	--truth-sensitivity-filter-level ${indel_filter_level} \
	--create-output-variant-index true \
	-mode INDEL
else
    cp $indels_recal $out_dir/tmp.indel.recalibrated.vcf
fi

recalibrated_vcf="$out_dir/joint_gvcfs.filtered.vcg.gz"
if [ -e $snps_recal -a ! -s $recalibrated_vcf ] ; then
   echo "Applying SNPs recalibration"
    $gatk_exec --java-options "-Xmx5g -Xms5g" \
	ApplyVQSR \
	-O ${recalibrated_vcf} \
	-V $out_dir/tmp.indel.recalibrated.vcf \
	--recal-file ${snps_recal} \
	--tranches-file ${snps_tranches} \
	--truth-sensitivity-filter-level ${snp_filter_level} \
	--create-output-variant-index true \
	-mode SNP
fi


# Collect variant calling metrics
metrics_prefix="$out_dir/joint_gvcfs"

    #java -Xmx6g -Xms6g -jar /usr/gitc/picard.jar \
    $gatk_exec --java-options "-Xmx6g -Xms6g" \
      CollectVariantCallingMetrics \
      --INPUT ${recalibrated_vcf} \
      --DBSNP ${dbSNP_vcf} \
      --SEQUENCE_DICTIONARY ${ref_dict} \
      --OUTPUT ${metrics_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${intervals}
