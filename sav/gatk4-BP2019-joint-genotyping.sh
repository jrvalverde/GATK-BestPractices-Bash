#!/bin/bash

################### THESE PROVIDE DEFAULT VALUES ###################
################## OVERRIDEN IN gatk4-BP2019.in.sh #################
data_dir="./hg19"
gatk_bundle_dir="$data_dir/gatk-bundle"
gatk_exec="$HOME/contrib/gatk4/gatk"
gatk_path="$HOME/contrib/gatk4/"
gotc_path="$HOME/contrib/gatk4/"
gitc_path="$HOME/contrib/gatk3/"
align_dir='align.gatk'
gvcf_dir='g.vcf'
workspace_dir_name='joint_gvcf'
out_dir='joint_gvcf_analysis'
#
ref_name=ucsc.hg19
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"
hapmap_vcf="$gatk_bundle_dir/hapmap_3.3.hg19.sites.vcf"
omni_vcf="$gatk_bundle_dir/1000G_omni2.5.hg19.sites.vcf"
Mills_vcf="$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
G1000_vcf="$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"
#
SNP_VQSR_downsampleFactor=10
snp_filter_level=99.7
indel_filter_level=99.7
indel_recalibration_annotation_values=("FS" "ReadPosRankSum" "MQRankSum" "QD" "SOR" "DP")
snp_recalibration_annotation_values=("QD" "MQRankSum" "ReadPosRankSum" "FS" "MQ" "SOR" "DP")
#MG=6
#MG=4
#MG=3
MG=2	# Max-gaussians parameter for VariantRecalibrator.
        # The --max-gaussians parameter sets the expected number of clusters
	# in modeling. If a dataset gives fewer distinct clusters, e.g. as can 
        # happen for smaller data, then the tool will tell you there is 
        # insufficient data with a No data found error message. In this case, 
        # try decrementing the --max-gaussians value. 

# get actual definitions from a local file
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




ref_fasta="$data_dir/ucsc.hg19.fasta"
ref_dict="$data_dir/ucsc.hg19.dict"

uBam_list=$flowcell_unmapped_bams_list
picard_interval_list=$study_interval_list
sample_name_map='sample_name.map'

function create_sample_name_map() {
    ### NOTE that these guards, all throughout this file,
    ### in all functions, ensure that all the parameters
    ### have been supplied and, hence, render needless any
    ### subsequent default assignment for missing parameters
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} uBam.list"
    else
        echo ">>> ${FUNCNAME[0]} $1 "
    fi
    uBam_list=$1
    sample_name_map='sample_name.map'

    #rm sample_name.map
    # prepare sample name map from uBam.list
    if [ ! -e $sample_name_map ] ; then
        echo ">>> Creating sample_name.map"
        # THIS SHOULD MATCH EXACTLY THE VALUES USED TO GENERATE
        # SAMPLE NAMES FOR THE ORIGINAL UNMAPPED BAM FILES
        # CHECK 00-convert-fastq-to-ubam.sh
        cat $uBam_list | while read bam ; do
            F=`basename $bam .bam`	# remove bam extension
            sample_id=`echo $F | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
            #SAMPLE_NAME=${sample_id}-`echo $F | cut -d'-' -f1`
            SAMPLE_NAME=${sample_id}`echo $F | cut -d'-' -f1`
            file_name=`ls $gvcf_dir/$F.*.ann.g.vcf`
            echo "$SAMPLE_NAME	$file_name"
        done > $sample_name_map
    fi
    
}


# This will import all gvcfs into a subdirectory and will create a tar file
# with the subdirectory contents
function import_gvcfs()
{
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sample_map reference interval_list work_dir"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi
    sample_name_map=$1
    ref_fasta=$2
    intervals=$3
    workspace_dir_name=$4
    batch_size=0	# read all samples at once
    
    if [ ! -s ${workspace_dir_name}.tar ] ; then
        echo ">>> Importing GVCFs"

        rm -rf ${workspace_dir_name}

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
}

# This will use the contents of the workspace containing all imported
# GVCFs to create a joint file and a combined selection of variants
# using the specified output file name.
# The output files will be
#	$out_dir/joint_gvcf.vcf.gz
#	$out_dir/joint_gvcf_combined.g.vcf
#
function genotype_gvcfs() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} work_dir reference intervals out_dir"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3$ $4"
    fi
    WORKSPACE=$1
    ref_fasta=$2
    intervals=$3
    joint_vcf=$4

    
    if [ ! -s $joint_vcf ] ; then
        echo ">>> Genotyping GVCFs"

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

    combined_vcf=${joint_vcf%.vcf*}.combined.g.vcf
    if [ ! -s $combined_vcf ] ; then
        echo ">>> Combining GVCFs"
        # generate a flat multisample GVCF file from the GenomicsDB
        $gatk_exec SelectVariants \
            -R $ref_fasta \
            -V gendb://$WORKSPACE \
            -O $combined_vcf

    fi
}


# do a hard filtering of all variants
# if no arguments are provided, then the defaults will be used:
#	$out_dir/joint_gvcf.vcf.gz
#	$out_dir/joint_gvcf.variant_filtered.vcf.gz
function filter_variants() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} [joint_vcf] [filtered_vcf]"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    # this means: use value of var $1 if set, otherwise, use the
    # value after :-
    joint_vcf=${1:-$out_dir/joint_gvcf.vcf.gz}
    # create default output file name by modifying the suffix
    default_filtered=${joing_vcf/%.vcf.gz/.variant_filtered.vcf.gz}

    variant_filtered_vcf=${2:-$default_filtered}
    #variant_filtered_vcf=${2:-$out_dir/joint_gvcf.variant_filtered.vcf.gz}
    
    excess_het_thresold=54.69

    # Hard-filter a large cohort callset on ExcessHet using VariantFiltration
    if [ ! -s ${variant_filtered_vcf} ] ; then
        echo ">>> Hard filtering"
        $gatk_exec --java-options "-Xmx3g -Xms3g" \
            VariantFiltration \
            --filter-expression "ExcessHet > 54.69" \
            --filter-name ExcessHet \
            -O ${variant_filtered_vcf} \
            -V ${joint_vcf}
    #        --filter-expression "ExcessHet > ${excess_het_threshold}" \
    fi
    ## LPZ
    ##FILTER=<ID=HighHaplotypeScore,Description="'HaplotypeScore > 13.0'">
    ##FILTER=<ID=HighQualZero,Description="'(MQ0 / DP) > 0.1'">
    ##FILTER=<ID=IndelFilter,Description="'QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0'">
    ##FILTER=<ID=LowDepth,Description="Applied to variants characterized in regions with 10X">
    ##FILTER=<ID=LowDepthAlt,Description="Applied to the variants with a number of reads with the ALT allele less below a threshold">
    ##FILTER=<ID=LowGQ,Description="'GQ < 20'">
    ##FILTER=<ID=LowMQRankSum,Description="'MQRankSum < -12.5'">
    ##FILTER=<ID=LowQual,Description="Low quality">
    ##FILTER=<ID=LowQualByDepth,Description="'QD < 2.0'">
    ##FILTER=<ID=LowRMSMappingQual,Description="'MQ < 30.0'">
    ##FILTER=<ID=LowReadDepth,Description="'DP < 16'">
    ##FILTER=<ID=LowReadPosRankSum,Description="'ReadPosRankSum < -8.0'">
    ##FILTER=<ID=PASS,Description="All filters passed">
    ##FILTER=<ID=StrandBias,Description="'FS > 60.0'">
}


function make_nocall_file() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} var_filtered var_filtered_nocall"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi

    variant_filtered_vcf=${1:-$out_dir/joint_gvcf.variant_filtered.vcf.gz}
    # set default output by changing the ending .vcf.gz to .nocal.vcf.gz
    default_nocall=${variant_filtered_vcf/%.vcf.gz/.nocall.vcf.gz}

    # use name provided or the default extended extension if none    
    variant_filtered_nocall_vcf=${2:-$default_nocall}
    #variant_filtered_nocall_vcf=${2:-$out_dir/joint_gvcf.variant_filtered.nocall.vcf.gz}

    # transform filtered genotypes to no-call
    if [ ! -s $variant_filtered_nocall_vcf ] ; then
        echo ">>> Making nocall file"
        $gatk_exec SelectVariants \
        -V $variant_filtered_vcf \
        --set-filtered-gt-to-nocall \
        -O $variant_filtered_nocall_vcf
    fi
}



function make_sites_only_vcf() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} var_filtered var_filtered_sites_only"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    variant_filtered_vcf=${1:-$out_dir/joint_gvcf.variant_filtered.vcf.gz}
    default_sites=${variant_filtered_vcf/%.vcf.gz/.sites_only.vcf.gz}

    sites_only_vcf=${2:-$default_sites}
    #sites_only_vcf=${2:-$out_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    
    # Create sites-only VCF with MakeSitesOnlyVcf
    if [ ! -s ${sites_only_vcf} ] ; then
        echo ">>> Making sites-only VCF"
        #java -Xmx3g -Xms3g -jar $gitc_path/picard.jar \
        $gatk_exec --java-options "-Xmx3g -Xms3g" \
            MakeSitesOnlyVcf \
            -INPUT ${variant_filtered_vcf} \
            -OUTPUT ${sites_only_vcf}
    fi
}


function recalibrate_indels() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf recal_indels_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi

    sites_only_vcf=${1:-$out_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    default_indels=${sites_only/%.vcf.gz/.indels.recal.vcf}

    indels_recal=${2:-$default_indels}
    #indels_recal="$out_dir/joint_gvfcs.indels.recal.vcf"
    
    indels_model_report=${indels_recal/%recal.vcf/model.report}
    #indels_model_report="$out_dir/joint_gvfcs.indels.model.report"
    indels_tranches=${indels_recal/%recal.vcf/tranches}
    #indels_tranches="$out_dir/joint_gvfcs.indels.tranches"
    indels_R=${indels_recal/recal.vcf/R}
    #indels_R="$out_dir/joint_gvfcs.indels.R"


    # Prepare arguments with annotations to use as reference for the recalibration
    # (see config file)
    an_arg=''
    for i in ${indel_recalibration_annotation_values[@]} ; do
        an_arg="$an_arg -an $i"
    done

    if [ ! -s $indels_recal -o ! -s $indels_tranches ] ; then
        echo ">>> Recalibrating INDELS"

        # we do it all at once in one single step
        # thus, we do not need to save the model, but we'll do anyway
        $gatk_exec --java-options "-Xmx24g -Xms24g" \
            VariantRecalibrator \
            -V ${sites_only_vcf} \
            -O ${indels_recal} \
            --tranches-file ${indels_tranches} \
            --output-model ${indels_model_report} \
            --trust-all-polymorphic \
            -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
            `# -an FS` \
            `# -an ReadPosRankSum` \
            `# -an QD` \
            `# -an SOR` \
            `# -an DP` \
            $an_arg \
            -mode INDEL \
            --max-gaussians $MG \
            -resource:mills,known=false,training=true,truth=true,prior=12 \
	        ${Mills_vcf} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2 \
	        ${dbSNP_vcf} \
	    -L $intervals \
            #--rscript-file ${indels_R}
    fi

    # repeat but using more reference databases
    # we make the new names substituting indels by indels.+
    indels_recal_plus=${indels_recal//indels/indels.+}
    #indels_recal_plus="$out_dir/joint_gvfcs.indels.+.recal.vcf"
    indels_model_report_plus=${indels_model_report//indels/indels.+}
    #indels_model_report_plus="$out_dir/joint_gvfcs.indels.+.model.report"
    indels_tranches_plus=${indels_tranches//indels/indels.+}
    #indels_tranches_plus="$out_dir/joint_gvfcs.indels.+.tranches"
    indels_plus_R=${indels_R//indels/indels.+}
    #indels_plus_R="$out_dir/joint_gvfcs.indels.+.R"
    
    if [ ! -s $indels_recal_plus -o ! -s $indels_tranches_plus ] ; then
        echo ">>> Recalibrating INDELS: USING ADDITIONAL REFERENCE DATA!"

        $gatk_exec --java-options "-Xmx3g -Xms3g" \
                VariantRecalibrator \
                -V ${sites_only_vcf} \
                -O ${indels_recal_plus} \
                --tranches-file ${indels_tranches_plus} \
                --output-model ${indels_model_report_plus} \
                --trust-all-polymorphic \
                -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
                `# -an FS -an QD -an SOR` \
                `# -an MQ -an GQ` \
                `# -an DP` \
                `# -an MQRankSum -an ReadPosRankSum` \
                $an_arg \
                -mode INDEL \
                -L $intervals \
                --max-gaussians $MG \
                -resource:hapmap,known=false,training=true,truth=true,prior=15 \
        	    ${hapmap_vcf} \
                -resource:mills,known=false,training=true,truth=true,prior=12 \
                      ${Mills_vcf} \
                -resource:omni,known=false,training=true,truth=true,prior=12 \
        	    ${omni_vcf} \
                -resource:1000G,known=false,training=true,truth=false,prior=10 \
        	    ${G1000_vcf} \
                -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        	    ${dbSNP_vcf} \
                #--rscript-file ${indels_plus_R}
    fi
}


function recalibrate_snps() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf recal_snps_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    sites_only_vcf=${1:-$out_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    default_snps=${sites_only_vcf/%.vcf.gz/.snps.recal.vcf}
    
    snps_recal=${2:-$default_snps}

    #  Calculate VQSLOD tranches for SNPs using VariantRecalibrator 
    snps_model_report=${snps_recal/%recal.vcf/model.report}
    snps_tranches=${snps_recal/%recal.vcf/tranches}
    snp_output_R=${snps_recal/%recal.vcf/R}

    downsampleFactor=$SNP_VQSR_downsampleFactor

    # Caution: here we use an external variable!!!
    # we should think this all better
    num_gvcfs=`wc -l $sample_name_map | cut -d' ' -f1`

    # build argument list with annotations to use in recalibration
    # (see config file)
    an_arg=''
    for i in ${snp_recalibration_annotation_values[@]} ; do
        an_arg="$an_arg -an $i"
    done
    
    if [ $num_gvcfs -gt 1000 ] ; then
        # Build a down-sampled model first and then recalibrate using this model
        if [ ! -s $snps_model_report -o ! -s $snps_tranches ] ; then
        echo ">>> Creating model for SNPs"

        $gatk_exec --java-options "-Xmx24g -Xms24g" \
            VariantRecalibrator \
            -V ${sites_only_vcf} \
            -O $snps_recal \
            --tranches-file $snps_tranches \
            --output-model ${snps_model_report} \
            --trust-all-polymorphic \
            -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
            `# -an FS -an QD -an SOR` \
            `# -an MQ -an GQ` \
            `# -an DP` \
            `# -an ReadPosRankSum -an MQRankSum` \
            $an_arg \
            -mode SNP \
            --sample-every-Nth-variant ${downsampleFactor} \
            --max-gaussians $MG \
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

            echo ">>> Applying model for SNPs"
            # Recalibrate SNPs with  the previous model
            $gatk_exec --java-options "-Xmx3g -Xms3g" \
              VariantRecalibrator \
              -V ${sites_only_nocall_vcf} \
              --input-model $snps_model_report \
              -O ${snps_recal} \
              --tranches-file ${snps_tranches} \
              --trust-all-polymorphic \
              -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
              `# -an FS` \
              `# -an ReadPosRankSum` \
              `# -an MQRankSum` \
              `# -an QD` \
              `# -an MQ` \
              `# -an SOR` \
              `# -an DP` \
              $an_arg \
              -mode SNP \
              --max-gaussians $MG \
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
              #--rscript-file ${snp_output_R}
        fi
    else
        if [ ! -e "$snps_recal" -o ! -e "$snps_tranches" ] ; then
	    echo ">>> Recalibrating SNPs directly"
	    # Recalibrating SNPs directly without using an intermediate model
            # we still generate it, but all the work is done in one single step
            $gatk_exec --java-options "-Xmx3g -Xms3g" \
                VariantRecalibrator \
                -V ${sites_only_vcf} \
                -O ${snps_recal} \
                --output-model $snps_model_report \
                --tranches-file ${snps_tranches} \
                --trust-all-polymorphic \
                -tranche "100.0" -tranche "99.95" -tranche "99.9" -tranche "99.8" -tranche "99.6" -tranche "99.5" -tranche "99.4" -tranche "99.3" -tranche "99.0" -tranche "98.0" -tranche "97.0" -tranche "90.0" \
                `# -an FS -an ReadPosRankSum -an QD -an SOR -an DP` \
                $an_arg \
                -mode SNP \
                -L $intervals \
                --max-gaussians $MG \
                -resource:hapmap,known=false,training=true,truth=true,prior=15 \
        	    ${hapmap_vcf} \
                -resource:mills,known=false,training=true,truth=true,prior=12 \
                      ${Mills_vcf} \
                -resource:omni,known=false,training=true,truth=true,prior=12 \
        	    ${omni_vcf} \
                -resource:1000G,known=false,training=true,truth=false,prior=10 \
        	    ${G1000_vcf} \
                -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        	    ${dbSNP_vcf} \
                #--rscript-file ${snp_output_R}

        fi
    #-an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
    # MQ has no data and so MQRankSum has zero variance
    fi

}



# Apply Recalibration
function add_indel_vqsr() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf indels_vcf indels_tranches indels_vqsr"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi

    sites_vcf=${1:-$out_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    default_vqsr=${sites_vcf%.vcf*}.indels.recal.vqsr.vcf
    
    indels_recal=${2:-$out_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vcf}
    default_tranches=${indels_recal%recal.vcf*}tranches

    indels_tranches=${3:-default_tranches}
    indels_vqsr=${4:-default_vqsr}

    # take the sites only vcf and filter the indel recalibration at the desired 
    # level generating a temporary file
    if [ ! -e $indels_vqsr ] ; then
        echo ">>> Applying INDEL VQSR recalibration"

        $gatk_exec --java-options "-Xmx5g -Xms5g" \
	    ApplyVQSR \
	    -V ${sites_vcf} \
	    --recal-file ${indels_recal} \
	    -O ${indels_vqsr} \
	    --tranches-file ${indels_tranches} \
	    --truth-sensitivity-filter-level ${indel_filter_level} \
	    --create-output-variant-index true \
	    -mode INDEL
    fi
}


function add_snps_vqsr() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_vcf snps_vcf snps_tranches snps_vqsr"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi

    sites_vcf=${1:-$out_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vqsr.vcf}
    default_vqsr=${sites_vcf%.vcf*}.snps.recal.vqsr.vcf

    snps_recal=${2:-$out_dir/joint_gvfc.variant_filtered.sites_only.snps.recal.vcf}
    default_tranches=${snps_recal/%recal.vcf/tranches}

    snps_tranches=${3:-default_tranches}
    snps_vqsr=${4:-default_vqsr}

    # Take as input a (possibly annotated) sites file and filter
    # the recalibrated snps at the level required.

    if [ -e $snps_recal -a -e $snps_tranches -a ! -e ${snps_vqsr} ] ; then
       echo "Applying SNPs VQSR recalibration"
        $gatk_exec --java-options "-Xmx5g -Xms5g" \
	    ApplyVQSR \
	    -V ${sites_vcf} \
	    --recal-file ${snps_recal} \
	    -O ${snps_vqsr} \
	    --tranches-file ${snps_tranches} \
	    --truth-sensitivity-filter-level ${snp_filter_level} \
	    --create-output-variant-index true \
	    -mode SNP
    fi
}


# collect variant calling metrics
function collect_metrics() {
    if [ $# -ne 3 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf intervals_list metrix_prefix"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3"
    fi

    recalibrated_vcf=${1:-$out_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vqsr.snps.recal.vqsr.vcf}
    default_metrics=${1%vcf*}metrics

    intervals=$2

    metrics_prefix=${3-default_metrics}
    #metrics_prefix="$out_dir/joint_gvcfs_metrics"
    
    if [ ! -e "$recalibrated_vcf" ] ; then
        echo ">>> ${FUNCNAME[0]}: $recalibrated_vcf doesn't exist"
        return 1
    fi

    if [ ! -s ${metrics_prefix}.variant_calling_detail_metrics ] ; then
        #java -Xmx6g -Xms6g -jar /usr/gitc/picard.jar \
        $gatk_exec --java-options "-Xmx6g -Xms6g" \
          CollectVariantCallingMetrics \
          --INPUT ${recalibrated_vcf} \
          --DBSNP ${dbSNP_vcf} \
          --SEQUENCE_DICTIONARY ${ref_dict} \
          --OUTPUT ${metrics_prefix} \
          --THREAD_COUNT 8 \
          --TARGET_INTERVALS ${intervals}
    fi
}


if [ ! -d $out_dir ] ; then
    mkdir -p $out_dir
fi
# Make names by successively removing the trailing .vcf or .vcf.gz
# and adding a new suffix. This will generate very long names, but
# will also allow us to track the progress of the calculations.
# ${var%vcf*} will remove the shortest match of 'vcf*' from the
# end of the string (e.g. kk.vcf , kk.vcf.gz or kk.vcf-something 
# will all yield kk.)

                                                  joint_gvcf=$out_dir/joint_gvcf.vcf.gz
                                         joint_gvcf_filtered=${joint_gvcf%vcf*}variant_filtered.vcf.gz
                              joint_gvcf_filtered_sites_only=${joint_gvcf_filtered%vcf*}sites_only.vcf.gz

                                  joint_gvcf_filtered_nocall=${joint_gvcf_filtered%vcf*}nocall.vcf.gz
                       joint_gvcf_filtered_nocall_sites_only=${joint_gvcf_filtered_nocall%vcf*}sites_only.vcf.gz

                       joint_gvcf_filtered_sites_only_indels=${joint_gvcf_filtered_sites_only%vcf*}indels.recal.vcf
              joint_gvcf_filtered_sites_only_indels_tranches=${joint_gvcf_filtered_sites_only%vcf*}indels.tranches

                         joint_gvcf_filtered_sites_only_snps=${joint_gvcf_filtered_sites_only%vcf*}snps.recal.vcf
                joint_gvcf_filtered_sites_only_snps_tranches=${joint_gvcf_filtered_sites_only%vcf*}snps.tranches

                  joint_gvcf_filtered_sites_only_indels_vqsr=${joint_gvcf_filtered_sites_only_indels%vcf*}vqsr.vcf
        joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr=${joint_gvcf_filtered_sites_only_indels_vqsr%vcf*}snps.recal.vqsr.vcf
joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr_metrics=${joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr%vcf*}metrics


create_sample_name_map $uBam_list

import_gvcfs $sample_name_map $ref_fasta $picard_interval_list \
		$workspace_dir_name

genotype_gvcfs $workspace_dir_name $ref_fasta $picard_interval_list \
 		$joint_gvcf

filter_variants $joint_gvcf \
		$joint_gvcf_filtered

# make sites only for the filtered file
make_sites_only_vcf $joint_gvcf_filtered \
		$joint_gvcf_filtered_sites_only

make_nocall_file $joint_gvcf_filtered \
		$joint_gvcf_filtered_nocall

# make sites only for the filtered nocall file
make_sites_only_vcf $joint_gvcf_filtered_nocall \
		$joint_gvcf_filtered_nocall_sites_only

recalibrate_indels $joint_gvcf_filtered_sites_only \
		$joint_gvcf_filtered_sites_only_indels

recalibrate_snps $joint_gvcf_filtered_sites_only \
		$joint_gvcf_filtered_sites_only_snps

add_indel_vqsr $joint_gvcf_filtered_sites_only \
		$joint_gvcf_filtered_sites_only_indels \
                $joint_gvcf_filtered_sites_only_indels_tranches \
                $joint_gvcf_filtered_sites_only_indels_vqsr

add_snps_vqsr $joint_gvcf_filtered_sites_only_indels_vqsr \
		$joint_gvcf_filtered_sites_only_snps \
                $joint_gvcf_filtered_sites_only_snps_tranches \
                $joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr

collect_metrics $joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr        	\
		$picard_interval_list                                        	\
                $joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr_metrics


###   ###   ###
# repeat for indels_plus
#
# When we have called recalibrate_indels above,
# we have generated actually two files using different information.
# It is possible that one is better than the other. We'll therefore do 
# the summary using either of them. This will create the corresponding 
# names substituting "indels" by "indels.+"
indels_plus=${joint_gvcf_filtered_sites_only_indels//indels/indels.+}
indels_plus_tranches=${joint_gvcf_filtered_sites_only_indels_tranches//indels/indels.+}
indels_plus_vqsr=${joint_gvcf_filtered_sites_only_indels_vqsr//indels/indels.+}
indels_plus_vqsr_snps_vqsr=${joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr//indels/indels.+}
indels_plus_vqsr_snps_vqsr_metrics=${joint_gvcf_filtered_sites_only_indels_vqsr_snps_vqsr_metrics//indels/indels.+}

add_indel_vqsr  $joint_gvcf_filtered_sites_only \
		$indels_plus \
                $indels_plus_tranches \
                $indels_plus_vqsr

add_snps_vqsr   $indels_plus_vqsr \
		$joint_gvcf_filtered_sites_only_snps \
                $joint_gvcf_filtered_sites_only_snps_tranches \
                $indels_plus_vqsr_snps_vqsr

collect_metrics $indels_plus_vqsr_snps_vqsr        	\
		$picard_interval_list                                        	\
                $indels_plus_vqsr_snps_vqsr_metrics

