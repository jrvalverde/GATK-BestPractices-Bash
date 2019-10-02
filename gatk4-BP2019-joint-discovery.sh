## Copyright Broad Institute, 2018
## 
## This WDL implements the joint discovery and VQSR filtering portion of the GATK 
## Best Practices (June 2016) for germline SNP and Indel discovery in human 
## whole-genome sequencing (WGS) and exome sequencing data.
##
## Requirements/expectations :
## - One or more GVCFs produced by HaplotypeCaller in GVCF mode 
## - Bare minimum 1 WGS sample or 30 Exome samples. Gene panels are not supported.
##
## Outputs :
## - A VCF file and its index, filtered using variant quality score recalibration 
##   (VQSR) with genotypes for all samples present in the input VCF. All sites that 
##   are present in the input VCF are retained; filtered sites are annotated as such 
##   in the FILTER field.
##
## Note about VQSR wiring :
## The SNP and INDEL models are built in parallel, but then the corresponding 
## recalibrations are applied in series. Because the INDEL model is generally ready 
## first (because there are fewer indels than SNPs) we set INDEL recalibration to 
## be applied first to the input VCF, while the SNP model is still being built. By 
## the time the SNP model is available, the indel-recalibrated file is available to 
## serve as input to apply the SNP recalibration. If we did it the other way around, 
## we would have to wait until the SNP recal file was available despite the INDEL 
## recal file being there already, then apply SNP recalibration, then apply INDEL 
## recalibration. This would lead to a longer wall clock time for complete workflow 
## execution. Wiring the INDEL recalibration to be applied first solves the problem.
##
## Cromwell version support 
## - Successfully tested on v31
## - Does not work on versions < v23 due to output syntax
##
## Runtime parameters are optimized for Broad's Google Cloud Platform implementation. 
## For program versions, see docker containers. 
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the docker 
## page at https://hub.docker.com/r/broadinstitute/genomes-in-the-cloud/ for detailed
## licensing information pertaining to the included programs.

# for debugging
CALL=''
#CALL=echo

function JointGenotyping() {
  if [ $# -gt 1 ] ; then
      echo "Usage: ${FUNCNAME[0]} [gvcf_dir]"
      return
  elif [ $# -eq 1 ] ; then
      gvcf_dir=$1
      echo ">>> ${FUNCNAME[0]} GVCF_DIR OVERRIDEN! Using $gvcf_dir"
  #else
      #gvcf_dir=$gvcf_dir 	# a predefined global variable
  fi
  echo ">>> ${FUNCNAME[0]} $*"
   
  # Reference and Resources
  #global File ref_fasta
  #global File ref_dict

  #global File dbSNP_vcf
  
  #global Array[String] snp_recalibration_tranche_values
  #global Array[String] snp_recalibration_annotation_values
  #global Array[String] indel_recalibration_tranche_values
  #global Array[String] indel_recalibration_annotation_values

  #global File eval_interval_list
  #global File hapmap_resource_vcf
  #global File omni_resource_vcf
  #global File one_thousand_genomes_resource_vcf
  #global File mills_resource_vcf
  #global File axiomPoly_resource_vcf
  #global File dbsnp_resource_vcf = dbsnp_vcf
  
  #global File unpadded_intervals_file

  # ExcessHet is a phred-scaled p-value. We want a cutoff of anything more extreme
  # than a z-score of -4.5 which is a p-value of 3.4e-06, which phred-scaled is 54.69
  #global Float excess_het_threshold = 54.69
  #global Float snp_filter_level
  #global Float indel_filter_level
  #global Int SNP_VQSR_downsampleFactor

  #Array[String] unpadded_intervals = read_lines(DynamicallyCombineIntervals.output_intervals)


    if [ ! -d $joint_gvcf_out_dir ] ; then
        mkdir -p $joint_gvcf_out_dir
    fi

    sample_name_map="sample_name.map"
    create_sample_name_map_from_GVCF_dir $gvcf_dir $sample_name_map

    unpadded_intervals=$bait_interval_list
    
    $CALL ImportGVCFs $sample_name_map $unpadded_intervals $joint_gvcf_dir
      #input:
      #  sample_names = sample_names,
      #  interval = unpadded_intervals[idx],
      #  workspace_dir_name = "genomicsdb",
      #  input_gvcfs = input_gvcfs,
      #  input_gvcfs_indices = input_gvcfs_indices
      #  batch_size = 50
    joint_gvcf_tar="${joint_gvcf_dir}.tar"


    if [ ! -d $joint_gvcf_out_dir ] ; then
    	mkdir -p $joint_gvcf_out_dir
    fi

    $CALL GenotypeGVCFs $joint_gvcf_dir $unpadded_intervals
      #input:
      #  workspace_tar = ImportGVCFs.output_genomicsdb,
      #  interval = unpadded_intervals[idx],
      #  output_vcf_filename = "output.vcf.gz",
      #  ref_fasta = ref_fasta,
      #  ref_fasta_index = ref_fasta_index,
      #  ref_dict = ref_dict,
      #  dbsnp_vcf = dbsnp_vcf,
      #  dbsnp_vcf_index = dbsnp_vcf_index
    joint_genotyped_vcf=$joint_gvcf_out_dir/joint_gvcf.genotyped.vcf.gz


    $CALL HardFilterAndMakeSitesOnlyVcf $joint_genotyed_vcf
      #input:
      #  vcf = GenotypeGVCFs.output_vcf,
      #  vcf_index = GenotypeGVCFs.output_vcf_index,
      #  excess_het_threshold = excess_het_threshold,
      #  variant_filtered_vcf_filename = callset_name + "." + idx + ".variant_filtered.vcf.gz",
      #  sites_only_vcf_filename = callset_name + "." + idx + ".sites_only.variant_filtered.vcf.gz",
    joint_genotyped_filtered_vcf=${joint_genotyed_vcf%.vcf*}.variant_filtered.vcf.gz
    joint_genotyped_filtered_sites_only_vcf=${joint_genotyped_filtered_vcf%.vcf*}.sites_only.vcf.gz


    $CALL IndelsVariantRecalibrator $joint_genotyped_filtered_sites_only_vcf
    #input:
      #sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      #sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      #recalibration_filename = callset_name + ".indels.recal",
      #tranches_filename = callset_name + ".indels.tranches",
      #recalibration_tranche_values = indel_recalibration_tranche_values,
      #recalibration_annotation_values = indel_recalibration_annotation_values,
      #mills_resource_vcf = mills_resource_vcf,
      #mills_resource_vcf_index = mills_resource_vcf_index,
      #axiomPoly_resource_vcf = axiomPoly_resource_vcf,
      #axiomPoly_resource_vcf_index = axiomPoly_resource_vcf_index,
      #dbsnp_resource_vcf = dbsnp_resource_vcf,
      #dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
    indels_recalibration=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.indels.recal
    indels_tranches=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.indels.tranches
    indels_rscript=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.indels.R

  num_gvcfs=`GetNumberOfSamples $sample_name_map`
  if [ num_gvcfs -gt 10000 ] ; then
      $CALL SNPsVariantRecalibratorCreateModel $joint_genotyped_filtered_sites_only_vcf
      #input:
      #  sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      #  sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      #  recalibration_filename = callset_name + ".snps.recal",
      #  tranches_filename = callset_name + ".snps.tranches",
      #  recalibration_tranche_values = snp_recalibration_tranche_values,
      #  recalibration_annotation_values = snp_recalibration_annotation_values,
      #  downsampleFactor = SNP_VQSR_downsampleFactor,
      #  model_report_filename = callset_name + ".snps.model.report",
      #  hapmap_resource_vcf = hapmap_resource_vcf,
      #  hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      #  omni_resource_vcf = omni_resource_vcf,
      #  omni_resource_vcf_index = omni_resource_vcf_index,
      #  one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      #  one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      #  dbsnp_resource_vcf = dbsnp_resource_vcf,
      #  dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
    snps_recalibration=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.recal
    snps_tranches=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.tranches
    snps_model=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.model
    snps_rscript=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.R
      

    $CALL SNPsVariantRecalibrator $joint_genotyped_filtered_sites_only_vcf
      #input:
      #  sites_only_variant_filtered_vcf = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf[idx],
      #  sites_only_variant_filtered_vcf_index = HardFilterAndMakeSitesOnlyVcf.sites_only_vcf_index[idx],
      #  recalibration_filename = callset_name + ".snps." + idx + ".recal",
      #  tranches_filename = callset_name + ".snps." + idx + ".tranches",
      #  recalibration_tranche_values = snp_recalibration_tranche_values,
      #  recalibration_annotation_values = snp_recalibration_annotation_values,
      #  model_report = SNPsVariantRecalibratorCreateModel.model_report,
      #  hapmap_resource_vcf = hapmap_resource_vcf,
      #  hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      #  omni_resource_vcf = omni_resource_vcf,
      #  omni_resource_vcf_index = omni_resource_vcf_index,
      #  one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      #  one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      #  dbsnp_resource_vcf = dbsnp_resource_vcf,
      #  dbsnp_resource_vcf_index = dbsnp_resource_vcf_index
    snps_recalibration=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.recal
    snps_tranches=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.tranches
    snps_model=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.model
    snps_rscript=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.R
      
      
  else  
  #if (num_gvcfs <= 10000){
    $CALL SNPsVariantRecalibrator $joint_genotyped_filtered_sites_only_vcf
      #input:
      #    sites_only_variant_filtered_vcf = SitesOnlyGatherVcf.output_vcf,
      #    sites_only_variant_filtered_vcf_index = SitesOnlyGatherVcf.output_vcf_index,
      #    recalibration_filename = callset_name + ".snps.recal",
      #    tranches_filename = callset_name + ".snps.tranches",
      #    recalibration_tranche_values = snp_recalibration_tranche_values,
      #    recalibration_annotation_values = snp_recalibration_annotation_values,
      #    hapmap_resource_vcf = hapmap_resource_vcf,
      #    hapmap_resource_vcf_index = hapmap_resource_vcf_index,
      #    omni_resource_vcf = omni_resource_vcf,
      #    omni_resource_vcf_index = omni_resource_vcf_index,
      #    one_thousand_genomes_resource_vcf = one_thousand_genomes_resource_vcf,
      #    one_thousand_genomes_resource_vcf_index = one_thousand_genomes_resource_vcf_index,
      #    dbsnp_resource_vcf = dbsnp_resource_vcf,
      #    dbsnp_resource_vcf_index = dbsnp_resource_vcf_index,
    snps_recalibration=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.recal
    snps_tranches=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.tranches
    snps_model=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.model
    snps_rscript=${joint_genotyped_filtered_sites_only_vcf%.vcf*}.snps.R


    $CALL ApplyRecalibration $joint_genotyped_filtered_sites_only_vcf
      #input:
      #  recalibrated_vcf_filename = callset_name + ".filtered." + idx + ".vcf.gz",
      #  input_vcf = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf[idx],
      #  input_vcf_index = HardFilterAndMakeSitesOnlyVcf.variant_filtered_vcf_index[idx],
      #  indels_recalibration = IndelsVariantRecalibrator.recalibration,
      #  indels_recalibration_index = IndelsVariantRecalibrator.recalibration_index,
      #  indels_tranches = IndelsVariantRecalibrator.tranches,
      #  snps_recalibration = if defined(SNPsVariantRecalibratorScattered.recalibration) then select_first([SNPsVariantRecalibratorScattered.recalibration])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration]),
      #  snps_recalibration_index = if defined(SNPsVariantRecalibratorScattered.recalibration_index) then select_first([SNPsVariantRecalibratorScattered.recalibration_index])[idx] else select_first([SNPsVariantRecalibratorClassic.recalibration_index]),
      #  snps_tranches = select_first([SNPGatherTranches.tranches, SNPsVariantRecalibratorClassic.tranches]),
      #  indel_filter_level = indel_filter_level,
      #  snp_filter_level = snp_filter_level
    joint_genotyped_filtered_sites_only_indels_recal=${joint_genotyped_filtered_sites_only%.vcf*}.indels.recalibrated.vcf.gz
    joint_genotyped_filtered_sites_only_indels_recal_snps_recal=${joint_genotyped_filtered_sites_only_indels_recal%.vcf*}.snps.recalibrated.vcf.gz
    

    $CALL CollectVariantCallingMetrics \
        $joint_genotyped_filtered_sites_only_indels_recal_snps_recal \
        $unpadded_intervals
      #input:
      #  input_vcf = FinalGatherVcf.output_vcf,
      #  input_vcf_index = FinalGatherVcf.output_vcf_index,
      #  metrics_filename_prefix = callset_name,
      #  dbsnp_vcf = dbsnp_vcf,
      #  dbsnp_vcf_index = dbsnp_vcf_index,
      #  interval_list = eval_interval_list,
      #  ref_dict = ref_dict


  #output {
  #  outputs from the small callset path through the wdl
  #  File? output_vcf = FinalGatherVcf.output_vcf
  #  File? output_vcf_index = FinalGatherVcf.output_vcf_index
  #
  #  select metrics from the small callset path and the large callset path
  #  File detail_metrics_file = select_first([CollectMetricsOnFullVcf.detail_metrics_file, GatherMetrics.detail_metrics_file])
  #  File summary_metrics_file = select_first([CollectMetricsOnFullVcf.summary_metrics_file, GatherMetrics.summary_metrics_file])
  #
  #  output the interval list generated/used by this run workflow
  #  File output_intervals = DynamicallyCombineIntervals.output_intervals
  #}
}



# Obtain sample ID from the file name
# This function MUST be the same used to create the original unmapped
# BAM file.
function sample_id_from_fname() {
    echo "$1" | sed -e 's/.*_S/S/g' -e 's/_L.*//g'
}


# Obtainread group from the file name
# This function MUST be the same used to create the original unmapped
# BAM file.
function read_group_from_fname() {
    $prefix=$study_name
    sample_id=`echo $1 | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
    echo "$prefix${sample_id}"
}


# Obtain sample name from the file name
# This function MUST be the same used to create the original unmapped
# BAM file.
function sample_name_from_fname() {
    sample_id=`echo $1 | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
    echo ${sample_id}`echo $F | cut -d'-' -f1`
}

# Obtain platform unit from the fastQ file
# This function MUST be the same used to create the original unmapped
# BAM file.
function platform_unit_from_fastq() {
    sample_id=`echo $1 | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
    # we should check if it is compressed
    platform_unit=`zcat $1 | head -1 | cut -d':' -f1 | tr -d '@'`
    echo ${sample_id}${platform_unit}
}

function create_sample_name_map_from_GVCF_dir() {
    ### NOTE that these guards, all throughout this file,
    ### in all functions, ensure that all the parameters
    ### have been supplied and, hence, render needless any
    ### subsequent default assignment for missing parameters
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} gvcf_directory sample_name_map"
    else
        echo ">>> ${FUNCNAME[0]} $*"
    fi
    dir=$1
    sample_name_map=$2

    # check if already done or if update is needed
    if [ ! -s "$sample_name_map" -o \
        "$dir" -nt "$sample_name_map" ] ; then
        # ensure it works whether the files are compressed or not
        for i in $dir/*vcf.gz $dir/*vcf ; do
            SN=`sample_name_from_fname $i`
            echo -e "${SN}\t$i" >> "$sample_name_map"
        done
    fi
}


# no longer used, we try to avoid usin a uBam_list
function create_sample_name_map_from_uBam_list_v2() {
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
    if [ ! -e "$sample_name_map" -o \
        "$uBam_list" -nt "$sample_name_map" ] ; then
        echo ">>> Creating sample_name.map"
        # THIS SHOULD MATCH EXACTLY THE VALUES USED TO GENERATE
        # SAMPLE NAMES FOR THE ORIGINAL UNMAPPED BAM FILES
        # CHECK 00-convert-fastq-to-ubam.sh
        cat $uBam_list | while read bam ; do
            F=`basename $bam .bam`	# remove bam extension
            SAMPLE_NAME=`sample_name_from_fname $F`
            # it could end with .g.vcf .vcf .vcf.gz .g.vcf.gz ...
            file_name=`ls $gvcf_dir/$F*vcf*`
            echo "$SAMPLE_NAME	$file_name"
        done > $sample_name_map
    fi
    
}


# no longer used, we try to avoid usin a uBam_list
function create_sample_name_map_from_uBam_list_v1() {
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
    if [ ! -e "$sample_name_map"  -o \
        "$uBam_list" -nt "$sample_name_map" ] ; then
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


# will print the number of samples on stdout
function GetNumberOfSamples() {
  # Can't print messages because we produce our output on stdout
  # so we send it to stderr
  if [ $# -ne 1 ] ; then
      echo "Usage: ${FUNCNAME[0]} sample_name_map" 1>&2
  else
      echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4" 1>&2
  fi
  #File sample_name_map
  sample_name_map=$1

  (
    wc -l ${sample_name_map} | awk '{print $1}'
  )
  #output {
  #  Int sample_count = read_int(stdout())
  #}
}

# Import GVCFs into a joint-gvcf directory and make a tar file of them
function ImportGVCFs() {
  if [ $# -ne 3 ] ; then
      echo "Usage: ${FUNCNAME[0]} sample_map interval_list work_dir"
  else
      echo ">>> ${FUNCNAME[0]} $*"
  fi
  #Array[String] sample_names
  #Array[File] input_gvcfs
  #Array[File] input_gvcfs_indices
  #String interval
  #String workspace_dir_name
  sample_names=$1
  interval=$3
  joint_gvcf_dir=$4
  
  batch_size=0		# read all samples at once (default in GATK-BP was 50)
  if [ "$batch_size" -gt "100" ] ; then
      consolidate_arg="--consolidate true"
  else
      consolidate_arg=""	# default is false
  fi

  if [ ! -s "${joint_gvcf_dir}.tar" -o \
      "$sample_names" -nt "${joint_gvcf_dir}.tar" ] ; then
    # The memory setting here is very important and must be several GB lower
    # than the total memory allocated to the VM because this tool uses
    # a significant amount of non-heap memory for native libraries.
    # Also, testing has shown that the multithreaded reader initialization
    # does not scale well beyond 5 threads, so don't increase beyond that.
    ${gatk_exec} \
        --java-options "-Xmx4g -Xms4g" \
        GenomicsDBImport \
        --genomicsdb-workspace-path ${joint_gvcf_dir} \
        --batch-size ${batch_size} \
        $consolidate_arg \
        -L ${interval} \
        --sample-name-map $sample_names \
        --reader-threads 5 \
        --interval-padding 500

    # this is only needed if we split the job over many subjobs in parallel
    tar -cf ${joint_gvcf_dir}.tar ${joint_gvcf_dir}

  fi
  #output {
  #  File output_genomicsdb = "${workspace_dir_name}.tar"
  #}
}


# extract joint-gvcf tar file and genotype the variants
function GenotypeGVCFs() {
  if [ $# -ne 2 ] ; then
      echo "Usage: ${FUNCNAME[0]} joint_gvcf_dir intervals"
  else
      echo ">>> ${FUNCNAME[0]} $*"
  fi
  #File workspace_tar
  #String interval
  #String output_vcf_filename
  #File ref_fasta
  #File ref_dict
  #File dbsnp_vcf
  joint_gvcf_dir=$1
  ref_fasta=$2
  interval=$3
  output_vcf_filename=${joint_gvcf_out_dir}/joint_gvcf.genotyped.vcf.gz

  joint_gvcf_tar=${joint_gvcf_dir}.tar

#### CHECK ACTUAL NAME OF OUTPUT GENERATED
  if [ ! -e "$output_vcf_filename" -o \
      "$joint_gvcf_tar" -nt "$output_vcf_filename" ] ; then
    set -e

    # this is only needed if we splitted the job over many subjobs in parallel
    # to consolidate the output of all subjobs into a single directory
    #	but letting it here will show on-screen that it runs and give feedback
    #   plus it allows us to check if an update is needed
    tar -xf ${joint_gvcf_tar}
    
    ${gatk_exec} \
        --java-options "-Xmx5g -Xms5g" \
        GenotypeGVCFs \
        -R ${ref_fasta} \
        -O ${output_vcf_filename} \
        -D ${dbsnp_vcf} \
        -G StandardAnnotation \
        --only-output-calls-starting-in-intervals \
        --use-new-qual-calculator \
        -V gendb://$joint_gvcf_dir \
        -L ${interval}
  fi
  #output {
  #  File output_vcf = "${output_vcf_filename}"
  #  File output_vcf_index = "${output_vcf_filename}.tbi"
  #}
}


function HardFilterAndMakeSitesOnlyVcf() {
  if [ $# -ne 1 ] ; then
      echo "Usage: ${FUNCNAME[0]} joint_vcf"
  else
      echo ">>> ${FUNCNAME[0]} $1"
  fi
  #File vcf
  #Float excess_het_threshold
  #String variant_filtered_vcf_filename
  #String sites_only_vcf_filename

  vcf=$1
  variant_filtered_vcf_filename=${vcf%.vcf*}.variant_filtered.vcf.gz
  sites_only_vcf_filename=${variant_filtered_vcf_filename%.vcf*}.sites_only.vcf.gz
  
  # ensure we have a default value in case it is undefined
  excess_het_threshold=${excess_het_threshold:-"54.69"}

  if [ ! -e "$variant_filtered_vcf_filename" -o \
       "${vcf}" -nt "$variant_filtered_vcf_filename" ] ; then
    ${gatk_exec} \
      --java-options "-Xmx3g -Xms3g" \
      VariantFiltration \
      --filter-expression "ExcessHet > ${excess_het_threshold}" \
      --filter-name ExcessHet \
      -O ${variant_filtered_vcf_filename} \
      -V ${vcf}
  fi
  
  if [ ! -e "$sites_only_vcf_filename" -o \
       "${variant_filtered_vcf_filename}" -nt "${sites_only_vcf_filename}" ] ; then
    ${gatk_exec} --java-options "-Xmx3g -Xms3g" \
      MakeSitesOnlyVcf \
      --INPUT ${variant_filtered_vcf_filename} \
      --OUTPUT ${sites_only_vcf_filename}
  fi
  #output {
  #  File variant_filtered_vcf = "${variant_filtered_vcf_filename}"
  #  File variant_filtered_vcf_index = "${variant_filtered_vcf_filename}.tbi"
  #  File sites_only_vcf = "${sites_only_vcf_filename}"
  #  File sites_only_vcf_index = "${sites_only_vcf_filename}.tbi"
  #}
}


# Create recalibrator for INDELs using the sites_only VCF
function IndelsVariantRecalibrator() {
  if [ $# -ne 1 ] ; then
      echo "Usage: ${FUNCNAME[0]} sites_only_vcf"
  else
      echo ">>> ${FUNCNAME[0]} $1"
  fi
  #String recalibration_filename
  #String tranches_filename
  #Array[String] recalibration_tranche_values
  #Array[String] recalibration_annotation_values
  #File sites_only_variant_filtered_vcf
  #File mills_resource_vcf
  #File axiomPoly_resource_vcf
  #File dbsnp_resource_vcf
  
  sites_only_variant_filtered_vcf=$1

  recalibration_filename=${sites_only_variant_filtered_vcf%.vcf*}.indels.recal
  tranches_filename=${sites_only_variant_filtered_vcf%.vcf*}.indels.tranches
  rscript_filename=${sites_only_variant_filtered_vcf%.vcf*}.indels.R
  mills_resource_vcf="$Mills_vcf"
  axiomPoly_resource_vcf="$axiomPoly_vcf"
  dbsnp_resource_vcf="$dbSNP_vcf"

  # prepare special arguments
  an_arg=''
  for i in ${indel_recalibration_annotation_values[@]} ; do
      an_arg="$an_arg -an $i"
  done
  
  default_recalibration_tranche_values=(
      100 99.95 99.9 99.8 99.7 99.6 99.5 99.4 99.3 99.0 98.0 97.0 05.0 90.0
      )
  recalibration_tranche_values=${recalibration_tranche_values:-default_recalibration_tranche_values}
  tranche_arg=""
  for i in ${recalibration_tranche_values[@]} ; do
      tranche_arg="$tranche_arg -tranche $i"
  done
  
  axiomp_arg=''
  if [ "$axiomPoly_vcf" != "" -a -e "$axiomPoly_vcf" ] ; then
      axiom_arg="--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_resource_vcf}"
  fi


  if [ ! -e "$tranches_filename" -o ! -e "$recalibration_filename" -o \
      "${sites_only_variant_filtered_vcf}" -nt "$recalibration_filename" -o \
      "${sites_only_variant_filtered_vcf}" -nt "$tranches_filename" ] ; then
    ${gatk_exec} \
      --java-options "-Xmx24g -Xms24g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      $tranche_arg \
      $an_arg \
      -mode INDEL \
      --max-gaussians $MG \
      --resource:mills,known=false,training=true,truth=true,prior=12 \
      	${mills_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2 \
      	${dbsnp_resource_vcf} \
      $axiom_arg \
      --rscript-file ${rscript_filename}
  fi
  #output {
  #  File recalibration = "${recalibration_filename}"
  #  File recalibration_index = "${recalibration_filename}.idx"
  #  File tranches = "${tranches_filename}"
  #}
}


# For SNPs, when we have many VCFs (> 1000), we first create a model and 
# then (later) will use it to make a recalibrator
function SNPsVariantRecalibratorCreateModel() {
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1"
    fi
  #String recalibration_filename
  #String tranches_filename
  #Int downsampleFactor
  #String model_report_filename
  #global Array[String] recalibration_tranche_values
  #global Array[String] recalibration_annotation_values
  #global File sites_only_variant_filtered_vcf

  #File hapmap_resource_vcf
  #File omni_resource_vcf
  #File one_thousand_genomes_resource_vcf
  #File dbsnp_resource_vcf
  sites_only_variant_filtered_vcf="$1"
  
  recalibration_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.recal
  tranches_filename==${sites_only_variant_filtered_vcf%.vcf*}.snps.tranches
  model_report_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.model
  rscript_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.R

  hapmap_resource_vcf="$hapmap_vcf"
  omni_resource_vcf="$omni_vcf"
  one_thousand_genomes_resource_vcf="$G1000_vcf"
  dbsnp_resource_vcf="$dbSNP_vcf"
  
  downsampleFactor=$SNP_VQSR_downsampleFactor

  an_arg=''
  for i in ${snp_recalibration_annotation_values[@]} ; do
      an_arg="$an_arg -an $i"
  done

  default_recalibration_tranche_values=(
      100 99.95 99.9 99.8 99.7 99.6 99.5 99.4 99.3 99.0 98.0 97.0 05.0 90.0
      )
  recalibration_tranche_values=${recalibration_tranche_values:-default_recalibration_tranche_values}
  tranche_arg=""
  for i in ${recalibration_tranche_values[@]} ; do
      tranche_arg="$tranche_arg -tranche $i"
  done


  if [ ! -e "$recalibration_filename" -o \
       ! -e "$tranches_filename" -o \
       ! -e "$model_report_filename" -o \
       "$sites_only_variant_filtered_vcf" -nt "$recalibration_filename" -o \
       "$sites_only_variant_filtered_vcf" -nt "$tranches_filename" -o \
       "$sites_only_variant_filtered_vcf" -nt "$model_report_filename" ] ; then
  
    ${gatk_exec} \
      --java-options "-Xmx100g -Xms100g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      $tranche_arg \
      $an_arg \
      -mode SNP \
      --sample-every-Nth-variant ${downsampleFactor} \
      --output-model ${model_report_filename} \
      --max-gaussians $MG \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 \
      		${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12 \
      		${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 \
      		${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 \
      		${dbsnp_resource_vcf} \
      --rscript-file ${rscript_filename}
  fi
  #output {
  #  File model_report = "${model_report_filename}"
  #}
}


function SNPsVariantRecalibrator() {
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1"
    fi
  #String recalibration_filename
  #String tranches_filename
  #File? model_report [optional]
  #File sites_only_variant_filtered_vcf

  #Array[String] recalibration_tranche_values
  #Array[String] recalibration_annotation_values
  #File hapmap_resource_vcf
  #File omni_resource_vcf
  #File one_thousand_genomes_resource_vcf
  #File dbsnp_resource_vcf

  sites_only_variant_filtered_vcf=$1
  recalibration_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.recal
  tranches_filename==${sites_only_variant_filtered_vcf%.vcf*}.snps.tranches
  model_report_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.model
  rscript_filename=${sites_only_variant_filtered_vcf%.vcf*}.snps.out.R


  an_arg=''
  for i in ${snp_recalibration_annotation_values[@]} ; do
      an_arg="$an_arg -an $i"
  done

  default_recalibration_tranche_values=(
      100 99.95 99.9 99.8 99.7 99.6 99.5 99.4 99.3 99.0 98.0 97.0 05.0 90.0
      )
  recalibration_tranche_values=${recalibration_tranche_values:-default_recalibration_tranche_values}
  tranche_arg=""
  for i in ${recalibration_tranche_values[@]} ; do
      tranche_arg="$tranche_arg -tranche $i"
  done
  
  if [ -e "$model_report_filename" ] ; then
      # if there is a pre-computed model, use it
      model_arg="--input-model $model_report_filename"
      # in this case, we may expect other files to also exist, back them up???
      if [ -e "$recalibration_filename" ] ; then
          cp "$recalibration_filename" "$recalibration_filename".pre
      fi
      if [ -e "$tranches_filename" ] ; then
          cp "$tranches_filename" "$tranches_filename".pre
      fi
      if [ -e "$rscript_filename" ] ; then
          cp "$rscript_filename" "$rscript_filename".pre
      fi
  else
      # otherwise, save the one we create
      model_arg="--output-model $model_report_filename"
  fi


  if [ ! -e "$recalibration_filename" -o \
       ! -e "$tranches_filename" -o \
       ! -e "$model_report_filename" -o \
       "$sites_only_variant_filtered_vcf" -nt "$recalibration_filename" -o \
       "$sites_only_variant_filtered_vcf" -nt "$tranches_filename"  -o \
       "$sites_only_variant_filtered_vcf" -nt "$model_report_filename" ] ;  then

    ${gatk_exec} \
      --java-options "-Xmx3g -Xms3g" \
      VariantRecalibrator \
      -V ${sites_only_variant_filtered_vcf} \
      -O ${recalibration_filename} \
      --tranches-file ${tranches_filename} \
      --trust-all-polymorphic \
      $tranche_arg \
      $an_arg \
      -mode SNP \
      $model_arg \
      --max-gaussians $MG \
      --resource:hapmap,known=false,training=true,truth=true,prior=15 \
      		${hapmap_resource_vcf} \
      --resource:omni,known=false,training=true,truth=true,prior=12  \
      		${omni_resource_vcf} \
      --resource:1000G,known=false,training=true,truth=false,prior=10 \
      		${one_thousand_genomes_resource_vcf} \
      --resource:dbsnp,known=true,training=false,truth=false,prior=7 \
      		${dbsnp_resource_vcf} \
      --rscript-file ${rscript_filename}
  fi
  #output {
  #  File recalibration = "${recalibration_filename}"
  #  File recalibration_index = "${recalibration_filename}.idx"
  #  File tranches = "${tranches_filename}"
  #}
}


function ApplyRecalibration() {
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} input_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1"
    fi
  #String recalibrated_vcf_filename
  #File input_vcf
  #File indels_recalibration
  #File indels_tranches
  #File snps_recalibration
  #File snps_tranches

  #global Float indel_filter_level
  #global Float snp_filter_level

  input_vcf=$1
  indel_recal_filename=${input_vcf%.vcf*}.indels.recalibrated.vcf.gz
  recalibrated_vcf_filename=${indel_recal_filename%.vcf*}.snps.recalibrated.vcf.gz

  indels_recalibration=${sites_only_variant_filtered_vcf%.vcf*}.indels.recal
  indels_tranches=${sites_only_variant_filtered_vcf%.vcf*}.indels.tranches
  snps_recalibration=${sites_only_variant_filtered_vcf%.vcf*}.snps.recal
  snps_tranches=${sites_only_variant_filtered_vcf%.vcf*}.snps.tranches


  if [ ! -e 
 
    # Apply indel recalibration
    ${gatk_exec} \
    	--java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${indel_recal_filename} \
      -V ${input_vcf} \
      --recal-file ${indels_recalibration} \
      --tranches-file ${indels_tranches} \
      --truth-sensitivity-filter-level ${indel_filter_level} \
      --create-output-variant-index true \
      -mode INDEL

    # Apply snp recalibration
    ${gatk_exec} \
      --java-options "-Xmx5g -Xms5g" \
      ApplyVQSR \
      -O ${recalibrated_vcf_filename} \
      -V ${indel_recal_filename} \
      --recal-file ${snps_recalibration} \
      --tranches-file ${snps_tranches} \
      --truth-sensitivity-filter-level ${snp_filter_level} \
      --create-output-variant-index true \
      -mode SNP
  fi
  #output {
  #  File recalibrated_vcf = "${recalibrated_vcf_filename}"
  #  File recalibrated_vcf_index = "${recalibrated_vcf_filename}.tbi"
  #}
}


function CollectVariantCallingMetrics() {
  if [ $# -ne 2 ] ; then
      echo "Usage: ${FUNCNAME[0]} input_vcf interval_list"
  else
      echo ">>> ${FUNCNAME[0]} $1 $2"
  fi
  #File input_vcf
  #String metrics_filename_prefix
  #global File dbsnp_vcf
  #File interval_list
  #global File ref_dict

  input_vcf=$1
  
  metrics_filename_prefix=${input_vcf%.vcf*}.metrics
  dbsnp_vcf=$dbSNP_vcf

  if [ ! -e "$metrics_filename_prefix".variant_calling_detail_metrics -o \
       "$input_vcf" -nt "$metrics_filename_prefix".variant_calling_detail_metrics ] ; then
    ${gatk_path} --java-options "-Xmx6g -Xms6g" \
      CollectVariantCallingMetrics \
      --INPUT ${input_vcf} \
      --DBSNP ${dbsnp_vcf} \
      --SEQUENCE_DICTIONARY ${ref_dict} \
      --OUTPUT ${metrics_filename_prefix} \
      --THREAD_COUNT 8 \
      --TARGET_INTERVALS ${interval_list}
  fi
  #output {
  #  File detail_metrics_file = "${metrics_filename_prefix}.variant_calling_detail_metrics"
  #  File summary_metrics_file = "${metrics_filename_prefix}.variant_calling_summary_metrics"
  #}
}


JointGenotyping $*

