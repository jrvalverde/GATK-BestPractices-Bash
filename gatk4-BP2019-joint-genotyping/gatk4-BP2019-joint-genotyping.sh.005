#!/bin/bash

#mydir=`dirname "$0"`
mydir=`dirname "${BASH_SOURCE[0]}"`


################### THESE PROVIDE DEFAULT VALUES ###################
################## OVERRIDEN IN gatk4-BP2019.in.sh #################
# DEFINE DEFAULT CONFIGURATION OPTIONS
# ------------------------------------
# These are provided just in case the user configuration file does not
# contain them

## SAMPLE NAME AND UNMAPPED BAMS
ref_name="Homo_sapiens_assembly38"
study_name="SAMPLE"
flowcell_unmapped_bams_list="./uBam.list"
unmapped_bam_suffix=".bam"

## PROGRAM-SPECIFIC PARAMETERS
SNP_VQSR_downsampleFactor=10
snp_filter_level=99.7
indel_filter_level=99.7

# We should only include here annotations that we want in the report
# `# comment` is a bash-ism to include comments in line.
# the DP annotation invoked by Coverage) should not be used when working 
# with exome datasets (gatk docs, args for VQSR)
# LPZ
#-an "FS" -an "ReadPosRankSum" -an "MQRankSum" -an "QD" -an "SOR" -an "DP" \
#        -resource:axiomPoly,known=false,training=true,truth=false,prior=10 \
#	    ${axiomPoly_vcf} \
indel_recalibration_annotation_values=("FS" "QD" "SOR" "MQ" "GQ" "ReadPosRankSum" "MQRankSum" "DP")

# LPZ
#-an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
snp_recalibration_annotation_values=("FS" "QD" "SOR" "MQ" "GQ" "MQRankSum" "ReadPosRankSum" "DP")

#MG=6	# Max-gaussians parameter for VariantRecalibrator.
#MG=5	# The --max-gaussians parameter sets the expected number of clusters
#MG=4	# in modeling. If a dataset gives fewer distinct clusters, e.g. as can 
#MG=3	# happen for smaller data, then the tool will tell you there is
MG=2    # insufficient data with a No data found error message. In this case, 
        # try decrementing the --max-gaussians value. 

## DATA DIRECTORIES
data_dir="./hg38"
gatk_bundle_dir="$data_dir/gatk_bundle"
broad_reference="$data_dir/broad-reference/v0"
panel_dir='./panel'
align_dir='align.gatk'
gvcf_dir='g.vcf'
workspace_dir_name='joint_gvcf'
joint_gvcf_dir="joint_gvcf_analysis.mg=$MG"

## PATHS
tools='usr/bin'
bwa_path="${HOME}/contrib/bwa/bwa.kit"
gatk3=${HOME}/contrib/gatk3
gatk4=${HOME}/contrib/gatk4
gotc_path="${HOME}/contrib/gatk4/"
gitc_path="${HOME}/contrib/gatk3/"
gatk_exec="${HOME}/contrib/gatk4/gatk"
picard_exec="${HOME}/contrib/gatk4/gatk"
python_exec="${HOME}/contrib/anaconda3/bin/python"
VerifyBamID_home="${gatk4}/contrib/VerifyBamID2/"
VerifyBamID_exec="${VerifyBamID_home}/bin/VerifyBamID"

## AUXILIARY FILES
#bait_intervals="./panel/panel5.hg38.bait.intervals"
#bait_bed="./panel/panel5.hg38.bait.bed"
bait_interval_list="./panel/panel5.hg38.bait.interval_list"
#
# This one needs to be in GATK intervals (list) format so 
# we can add unmapped at the end
study_intervals="./panel/panel5.hg38.120.intervals"
#study_bed="./panel/panel5.hg38.120.bed"
# This may be a Picard interval_list, a GATK interval or list, or a BED file
# but a Picard interval_list is recommended
study_interval_list="./panel/panel5.hg38.120.interval_list"

## Broad Insitute contamination reference files
#  retrieved with gsutil -m cp -R gs://broad-references/hg38/v0. . 
contamination_sites_ud="${broad_reference}/Homo_sapiens_assembly38.contam.UD"
contamination_sites_bed="${broad_reference}/Homo_sapiens_assembly38.contam.bed"
contamination_sites_mu="${broad_reference}/Homo_sapiens_assembly38.contam.mu"
contamination_sites_V="${broad_reference}/Homo_sapiens_assembly38.contam.V"

## "KNOWN SITES RESOURCES"
dbSNP_vcf=( 
#	"$gatk_bundle_dir/dbsnp_138.hg38.vcf.gz"
#        "$gatk_bundle_dir/dbsnp_144.hg38.vcf.gz"
#        "$gatk_bundle_dir/dbsnp_146.hg38.vcf.gz"
	"$broad_reference/Homo_sapiens_assembly38.dbsnp138.vcf"
)
dbSNP_vcf="$broad_reference/Homo_sapiens_assembly38.dbsnp138.vcf"
G1000_vcf="$gatk_bundle_dir/1000G_phase1.snps.high_confidence.hg38.vcf.gz"
omni_vcf="$gatk_bundle_dir/1000G_omni2.5.hg38.vcf"
Mills_vcf="$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
axiomPoly_vcf="$gatk_bundle_dir/Axiom_Exome_Plus.genotypes.all_populations.poly.hg38.vcf.gz"
hapmap_vcf=( 
	"$gatk_bundle_dir/hapmap_3.3.hg38.vcf"
	"$gatk_bundle_dir/hapmap_3.3.grch38_pop_stratified_af.vcf"
)
hapmap_vcf="$gatk_bundle_dir/hapmap_3.3.hg38.vcf.gz"
known_indels_sites_VCFs=(
#    "$gatk_bundle_dir/dbsnp_138.hg38.vcf.gz"
#    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz"
    "${broad_reference}/Homo_sapiens_assembly38.known_indels.vcf.gz"
)


# LOAD USER CONFIGURATION OPTIONS
# -------------------------------
# get actual definitions from a local file IN THE WORKING DIRECTORY
# 
# I.e. this implies we need to have a config file in the directory from which
# we call this script.
#
if [ -s gatk4-BP2019.in.sh ] ; then
    . gatk4-BP2019.in.sh
else
    cp $my_dir/gatk4-BP2019.in.sh \
       ./gatk4-BP2019.in.sh
    echo "I have created a preferences file in this directory:"
    echo ""
    echo "        gatk4-BP2019.in.sh"
    echo ""
    echo "Please, open it in a text editor, adapt it to your needs and"
    echo "run this command again."
    exit
fi



# fill in other derived variables

ref_fasta="$data_dir/$ref_name.fasta"
ref_dict="$data_dir/$ref_name.dict"

uBam_list=$flowcell_unmapped_bams_list	# currently unused
picard_interval_list=$study_interval_list
sample_name_map='sample_name.map'


## sample_id_from_file()
#
# Usage:
#	sample_id=`sample_id_from_file $filename`
#
# Obtain the ID of the sample from the file it is stored in
#
#	This assumes that all the data corresponding to a single sample 
# is on a separate file, and that the file name contains different bits
# of information separated in identifiable "fields", among them, the
# information we want to retrieve.
#
#	For Sample-ID, we assume that fields are separated by '_' and 
# that the information is in a field tagged 'S', followed by another
# field tagged 'L'. Hence, we remove everything up to this field,
# identified by 'delimiter+label' (_S) leave the 'S' in, and then  
# remove everything from the next field (_L) to the end.
#
#	YMMV and you should adapt this function to provide a suitable
# sample identifier, possibly by lokking inside the file if needed.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function sample_id_from_file() {
    local F=`basename $1`
    
    echo "$F" | sed -e 's/.*_S/S/g' -e 's/_L.*//g'
}


## read_group_from_file
#
# Usage:
#	read_group=`read_group_from_file $filename`
#
# Obtain the read group for ALL reads in a file from the file it is stored in
#
#	This assumes that all the data corresponding to a single sample 
# is on a separate file, and that the file name contains different bits
# of information separated in identifiable "fields", among them, the
# information we want to retrieve.
#
#	For read-group, we assume that fields are separated by '_' and 
# that the information is in a field tagged 'S', followed by another
# field tagged 'L'. Hence, we remove everything up to this field,
# identified by 'delimiter+label' (_S) and then from the next field (_L)
# to the end. We combine this sample-id with the study name to define
# the read-group.
#
#	In other words, we use the sample-ID as read-group, but combine
# it with the study name to distinguish RG from SID.
#
#	YMMV and you should adapt this function to provide a suitable
# sample identifier, possibly by lokking inside the file if needed.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function read_group_from_file() {
    local F=`basename $1`
    local $prefix=$study_name		# 'study_name' must be globally defined
    local sample_id=`echo $F | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`

    echo "$prefix${sample_id}"
}


## sample_name_from_file
#
# Usage:
#	read_group=`sample_name_from_file $filename`
#
# Obtain the name of the sample from the file it is stored in
#
#	This assumes that all the data corresponding to a single sample 
# is on a separate file, and that the file name contains different bits
# of information separated in identifiable "fields", among them, the
# information we want to retrieve.
#
#	For read-group, we assume that fields are separated by '_' and 
# that the information is in a field tagged 'S', followed by another
# field tagged 'L'. Hence, we remove everything up to this field,
# identified by 'delimiter+label' (_S) and then from the next field (_L)
# to the end. We combine this sample-id with the first "field" in the
# file name to obtain the sample name. This first field is usually a unique
# identifier provided by the sequencer.
#
#	In other words, we use the sample-ID as read-group, but combine
# it with the sequencer sample ID to distinguish RG from SID.
#
#	YMMV and you should adapt this function to provide a suitable
# sample identifier, possibly by lokking inside the file if needed.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function sample_name_from_file() {
    local F=`basename $1`
    local sample_id=`echo $F | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`

    echo ${sample_id}`echo $F | cut -d'-' -f1`
}


## platform_unit_from_fastq
#
# Usage:
#	platform_unit=`platform_unit_from_fastq $fastq_filename`
#
# Obtain the identifier of the sequencing platform unit used from a
# fastq file.
#
#	We take the first line of the FastQ file and try to retrieve
# the platform unit identifier from the FastQ read header.
#
#	For extra safety WE combine the platform-unit with the sample
# id, so that this identifier may provide extra information. You may
# not want to do this. Also, we assume a specific organization of the
# read fastq header, you should check it applies in your own datasets.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function platform_unit_from_fastq() {
    local F=`basename $1`
    local sample_id=`echo $1 | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`
    # we should check if it is compressed
    local platform_unit=`zcat $1 | head -1 | cut -d':' -f1 | tr -d '@'`

    echo ${sample_id}${platform_unit}
}


## create_sample_name_map_from_GVCF_dir()
#
# Usage:
#	create_sample_name_map_from_GVCF_dir $gvcf_directory_or_folder
#
# Create a mapping from sample-names to corresponding GVCF files
#
#	We assume here that all GVCF files are gathered together in a
# single directory *that only contains* the GVCF files to analyze. If
# you have your GVCF files in a folder with more data files, it should
# be trivial to create a new directory with only the GVCF files you
# want. You may use links to reduce disk usage, e.g.
#
#	ls results
#	    file1.recalibrated.g.vcf
#	    file1.recalibrated.bam
#	    file1.vcf
#	    file2.recalibrated.g.vcf
#	    ...
#	mkdir g.vcf
#	cd g.vcf
#	ln -s ../results/*recalibrated.g.vcf .
#
#	This is so because joint genotyping uses a joint gvcf database
# that cannot be added to. It will only contain info on the gvcfs
# originally used to create it. If you later want to include a new
# gvcf in the analysis you must rebuild the database with all the
# old gvcfs plus the new one.
#
#	If you build the database indicating directly the GVCF names, then
# it will be difficult later on to know which files you used. If you first
# store the gvcfs in a specific directory and make the database, then you
# will know to which gvcf files the database corresponds to.
#
#	Later on, if you want to include new gvcfs in the analysis, all you
# need to do is add them to the gvcf directory/folder, remove the database
# and repeat the complete analysis (yes, you must repeat the complete
# analysis, and yes, I know folder and directory are the same, but many
# newbies ignore such a simple fact).
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function create_sample_name_map_from_GVCF_dir() {
    ### NOTE that these guards, all throughout this file,
    ### in all functions, ensure that all the parameters
    ### have been supplied and, hence, render needless any
    ### subsequent default assignment for missing parameters
    if [ $# -ne 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} gvcf_directory"
    else
        echo ">>> ${FUNCNAME[0]} $1 "
    fi
    local dir=$1
    local SN=''
    # global File sample_name_map

    # provide a default value, just in case
    sample_name_map=${sample_name_map:-'sample_name.map'}

    if [ ! -s "$sample_name_map" ] ; then
        # ensure it works whether the files are compressed or not
        for i in $dir/*vcf.gz $dir/*vcf ; do
            SN=`sample_name_from_file $i`
            echo -e "${SN}\t$i" >> "$sample_name_map"
        done
    fi
}


## create_sample_name_map_from_uBam_list_v2()
#
# Usage:
#	create_sample_name_map_from_uBam_list_v2 uBam.list
#
# Create a mapping from sample-names to gvcf files based on a list of
# original unmapped BAM files.
#
#	This function is intended to be used when the analysis is based
# on a file containing a list of unmapped BAMs (as in the standard Broad
# institute protocols). You do not need this in these scripts because we
# can take the file names (or a folder containing the files) form the
# command line, which is more convenient. But, should you prefer to work
# from a uBam.list file, this function will be handy to create the mapping
# required for joint genotyping (or joint variant discovery).
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
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
    local uBam_list=$1
    # global File sample_name_map
    
    # Provide a default value, just in case, and make it global
    sample_name_map=${sample_name_map:-'sample_name.map'}

    #rm sample_name.map
    # prepare sample name map from uBam.list
    if [ ! -e $sample_name_map ] ; then
        echo ">>> Creating sample_name.map"
        # THIS SHOULD MATCH EXACTLY THE VALUES USED TO GENERATE
        # SAMPLE NAMES FOR THE ORIGINAL UNMAPPED BAM FILES
        # CHECK 00-convert-fastq-to-ubam.sh
        cat $uBam_list | while read bam ; do
            F=`basename $bam .bam`	# remove bam extension
            SAMPLE_NAME=`sample_name_from_file $F`
            # it could end with .g.vcf .vcf .vcf.gz .g.vcf.gz ...
            file_name=`ls $gvcf_dir/$F*vcf*`
            echo "$SAMPLE_NAME	$file_name"
        done > $sample_name_map
    fi
    
}


## create_sample_name_map_from_uBam_list_v1()
#
# Usage:
#	create_sample_name_map_from_uBam_list_v1 uBam.list
#
# Create a mapping from sample-names to gvcf files based on a list of
# original unmapped BAM files.
#
#	This function is intended to be used when the analysis is based
# on a file containing a list of unmapped BAMs (as in the standard Broad
# institute protocols). You do not need this in these scripts because we
# can take the file names (or a folder containing the files) form the
# command line, which is more convenient. But, should you prefer to work
# from a uBam.list file, this function will be handy to create the mapping
# required for joint genotyping (or joint variant discovery).
#
#	What's the difference? "Use the source, Luke"
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
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

## import_gvcfs()
#
# Usage:
#	import_gvcfs sample_name_map reference interval_list work_dir
#
# This will import all gvcfs listed in 'sample_name_map' into a subdirectory 
# (work_dir) and will create a tar file with the subdirectory contents, using
# as reference fasta file 'reference', and considering only the sequence
# intervals listed in 'intrval_list' (padded with 500 nt on each side).
#
#	Currently we read all bam files at once. This may be excessive if
# you have too many files. In that case you may want to modify variable
# batch_size to a more suitable value like 50 (as in Broad's scripts).
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function import_gvcfs()
{
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sample_map reference interval_list work_dir"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi
    local sample_name_map=$1
    local ref_fasta=$2
    local intervals=$3
    local workspace_dir_name=$4
    #local batch_size=50	# read all samples at once <<<<<<< !!!!!!
    local batch_size=0	# read all samples at once <<<<<<< !!!!!!
    
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

	# NOTE: 
        # In the Broad scripts, this step is "parallelized" by splitting
        # the input in subsets, each subset generates a subdirectory,
        # and the subdirectory is archived with tar so that it may be
        # returned as a result; later on, all tar files will be extracted
        # to reconstruct the "whole" database working directory. 
        #
        # We do not split, so the tar file is not needed and indeed is
        # an inefficient step, but... IFF the tar file is not created,
        # then the execution condition should be changed to detect the
        # presence of the directory instead.
        #
        # NOTE:
        # Given the way the database directory is constructed in the Broad
        # scripts, one is left to ponder if the same trick could not be used
        # to add new samples to an analysis instead of rebuilding the whole
        # database. Maybe one of these days I'll have time to test it. Or
        # maybe the Broad's guys will realize themselves. Or maybe they 
        # will update the documentation to reflect it. Who knows?
        
    fi
}


## genotype_gvcfs()
#
# Usage:
#	genotype_gvcfs joint_db_dir ref_fasta intervals joint_gvcf
#
# This will use the contents of the workspace containing all imported
# GVCFs to create a joint file and a combined selection of variants
# using the specified output file name.
# The output files will be
#	$joint_gvcf_dir/joint_gvcf.vcf.gz
#	$joint_gvcf_dir/joint_gvcf_combined.g.vcf
#
# If the function is not used properly, it will output a message indicating
# the proper usage.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function genotype_gvcfs() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} work_dir reference intervals joint_gvcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3$ $4"
    fi
    WORKSPACE=$1
    ref_fasta=$2
    intervals=$3
    joint_vcf=$4	# the name of the output file

    if [ ! -d $WORKSPACE ] ; then
        if [ -e $WORKSPACE.tar ] ; then
            tar -xvf $WORKSPACE.tar
        elif [ -e $WORKSPACE.tar.gz ] ; then
            tar -xvf $WORKSPACE.tar.gz
        else
            echo ">>> ${FUNCNAME[0]} NO JOINT DATABASE $joint_vcf FOUND!"
            exit
        fi
    fi
    
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


## filter_variants()
#
# Usage:
#	filter_variants [joint_vcf] [filtered_vcf]
#
# do a hard filtering of all variants
#
#	Use the "joint_vcf" file provided, hard-filter it, and
# save the result as "filtered_vcf"
#
# 	If no arguments are provided, then these defaults will be used:
#		$joint_gvcf_dir/joint_gvcf.vcf.gz
#		$joint_gvcf_dir/joint_gvcf.variant_filtered.vcf.gz
#
#	Additional potential filters are listed as comments inside. 
# We will eventually transform this function to use a global array
# of filter-name / filter-condition pairs and build the arguments on 
# the fly, so that users can specify as many filters as they like in
# the global configuration file.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function filter_variants() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} [joint_vcf] [filtered_vcf]"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    # this means: use value of var $1 if set, otherwise, use the
    # value after :-
    local joint_vcf=${1:-$joint_gvcf_dir/joint_gvcf.vcf.gz}
    # create default output file name by modifying the suffix
    local default_filtered=${joing_vcf/%.vcf.gz/.variant_filtered.vcf.gz}

    local variant_filtered_vcf=${2:-$default_filtered}
    #variant_filtered_vcf=${2:-$joint_gvcf_dir/joint_gvcf.variant_filtered.vcf.gz}
    
    excess_het_thresold=${excess_het_threshold:-54.69}

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


## make_nocall_file
#
# Usage:
#	make_nocall_file filtered_vcf resulting_nocall_vcf
#
#	The nocall file will not be further used in the analysis 
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function make_nocall_file() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} var_filtered var_filtered_nocall"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi

    local variant_filtered_vcf=${1:-$joint_gvcf_dir/joint_gvcf.variant_filtered.vcf.gz}
    # set default output by changing the ending .vcf.gz to .nocal.vcf.gz
    local default_nocall=${variant_filtered_vcf/%.vcf.gz/.nocall.vcf.gz}

    # use name provided or the default extended extension if none    
    local variant_filtered_nocall_vcf=${2:-$default_nocall}
    #variant_filtered_nocall_vcf=${2:-$joint_gvcf_dir/joint_gvcf.variant_filtered.nocall.vcf.gz}

    # transform filtered genotypes to no-call
    if [ ! -s $variant_filtered_nocall_vcf ] ; then
        echo ">>> Making nocall file"
        $gatk_exec SelectVariants \
        -V $variant_filtered_vcf \
        --set-filtered-gt-to-nocall \
        -O $variant_filtered_nocall_vcf
    fi
}


## make_sites_only_vcf()
#
# Usage:
#	make_sites_only_vcf var_filtered var_filtered_sites_only
#
# Creates a VCF that contains all the site-level information for all records 
# in the input VCF but no genotype information.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function make_sites_only_vcf() {
    # Creates a VCF that contains all the site-level information for all records 
    # in the input VCF but no genotype information.
    #	=> some annotations (like GQ) will disappear
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} var_filtered var_filtered_sites_only"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    local variant_filtered_vcf=${1:-$joint_gvcf_dir/joint_gvcf.variant_filtered.vcf.gz}
    local default_sites=${variant_filtered_vcf/%.vcf.gz/.sites_only.vcf.gz}

    local sites_only_vcf=${2:-$default_sites}
    #sites_only_vcf=${2:-$joint_gvcf_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    
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


## recalibrate_indels()
#
# Usage:
#	recalibrate_indels sites_only_vcf recal_indels_vcf
#
# Build the INDEL recalibration model
#
# 	Note that tranches are hard-coded and the R-script file is not
# being currently generated. A newer, more general version of this
# function is available in script 'joint-discovery.sh'.
#
#	We make heavy use of global variables to specify the reference
# databases.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function recalibrate_indels() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf recal_indels_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi

    local sites_only_vcf=${1:-$joint_gvcf_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    local default_indels=${sites_only/%.vcf.gz/.indels.recal.vcf}

    local indels_recal=${2:-$default_indels}
    #indels_recal="$joint_gvcf_dir/joint_gvfcs.indels.recal.vcf"
    
    local indels_model_report=${indels_recal/%recal.vcf/model.report}
    #indels_model_report="$joint_gvcf_dir/joint_gvfcs.indels.model.report"
    local indels_tranches=${indels_recal/%recal.vcf/tranches}
    #indels_tranches="$joint_gvcf_dir/joint_gvfcs.indels.tranches"
    local indels_R=${indels_recal/recal.vcf/R}
    #indels_R="$joint_gvcf_dir/joint_gvfcs.indels.R"


    # Prepare arguments with annotations to use as reference for the recalibration
    # (see config file)
    local an_arg=''
    for i in ${indel_recalibration_annotation_values[@]} ; do
        an_arg="$an_arg -an $i"
    done
    local axiomp_arg=''
    if [ "$axiomPoly_vcf" != "" -a -e "$axiomPoly_vcf" ] ; then
        axiom_arg="--resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomPoly_vcf}"
    fi

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
            $axiom_arg \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2 \
	        ${dbSNP_vcf} \
	    -L $intervals \
            #--rscript-file ${indels_R}
    fi

    # UNDOCUMENTED FEATURE!!!
    #
    # repeat but using more reference databases
    #
    # we make the new names substituting indels by indels.+
    local indels_recal_plus=${indels_recal//indels/indels.+}
    #indels_recal_plus="$joint_gvcf_dir/joint_gvfcs.indels.+.recal.vcf"
    local indels_model_report_plus=${indels_model_report//indels/indels.+}
    #indels_model_report_plus="$joint_gvcf_dir/joint_gvfcs.indels.+.model.report"
    local indels_tranches_plus=${indels_tranches//indels/indels.+}
    #indels_tranches_plus="$joint_gvcf_dir/joint_gvfcs.indels.+.tranches"
    local indels_plus_R=${indels_R//indels/indels.+}
    #indels_plus_R="$joint_gvcf_dir/joint_gvfcs.indels.+.R"
    
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
    if [ ! -s $indels_recal_plus -o ! -s $indels_tranches_plus ] ; then
        echo ">>> WARNING: extra recalibration failed!"
    fi
}


## recalibrate_snps()
#
# Usage:
#	recalibrate_snps sites_only_vcf recal_snps_vcf
#
# Build recalibration model for SNPs
#
#	This function includes a double logic, one is for the case when
# the number of files to process is too large (> 1000), doing it in two
# separate steps, and the other is doing it all at once.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function recalibrate_snps() {
    if [ $# -ne 2 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf recal_snps_vcf"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2"
    fi
    local sites_only_vcf=${1:-$joint_gvcf_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    local default_snps=${sites_only_vcf/%.vcf.gz/.snps.recal.vcf}
    
    local snps_recal=${2:-$default_snps}

    #  Calculate VQSLOD tranches for SNPs using VariantRecalibrator 
    local snps_model_report=${snps_recal/%recal.vcf/model.report}
    local snps_tranches=${snps_recal/%recal.vcf/tranches}
    local snp_output_R=${snps_recal/%recal.vcf/R}

    local downsampleFactor=$SNP_VQSR_downsampleFactor

    # Caution: here we use an external variable!!!
    # we should think this all better
    local num_gvcfs=`wc -l $sample_name_map | cut -d' ' -f1`

    # build argument list with annotations to use in recalibration
    # (see config file)
    local an_arg=''
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
            -resource:omni,known=false,training=true,truth=false,prior=12.0 \
                $omni_vcf \
            -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
                $G1000_vcf \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
                $dbSNP_vcf \
	    -L $intervals \
            #-resource:mills,known=false,training=true,truth=true,prior=12.0 \
            #    ${Mills_vcf} \

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
              -resource:omni,known=false,training=true,truth=false,prior=12.0 \
                  $omni_vcf \
              -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
                  $G1000_vcf \
              -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
                  $dbSNP_vcf \
              #-resource:mills,known=false,training=true,truth=true,prior=12.0 \
              #    ${Mills_vcf} \
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
                -resource:omni,known=false,training=true,truth=true,prior=12 \
        	    ${omni_vcf} \
                -resource:1000G,known=false,training=true,truth=false,prior=10 \
        	    ${G1000_vcf} \
                -resource:dbsnp,known=true,training=false,truth=false,prior=7 \
        	    ${dbSNP_vcf} \
                #-resource:mills,known=false,training=true,truth=true,prior=12 \
                #      ${Mills_vcf} \
                #--rscript-file ${snp_output_R}

        fi
    #
    #-an FS -an ReadPosRankSum -an MQRankSum -an QD -an MQ -an SOR -an DP \
    #
    # NOTE:
    # MQ has no data and so MQRankSum has zero variance
    # Hence, MQRankSum should not be included in the recalibration annotation
    
    fi

}



## add_indel_vqsr()
#
# Usage: 
#	add_indel_vqsr sites_only_vcf indels_vcf indels_tranches indels_vqsr
#
# Apply INDEL Recalibration
#
# Use the INDEL recalibration model to add the INDEL VQSR to the specified
# VCF file
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function add_indel_vqsr() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf indels_vcf indels_tranches indels_vqsr"
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi

    local sites_vcf=${1:-$joint_gvcf_dir/joint_gvcf.variant_filtered.sites_only.vcf.gz}
    local indels_recal=${2:-$joint_gvcf_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vcf}

    local default_vqsr=${sites_vcf%.vcf*}.indels.recal.vqsr.vcf
    local default_tranches=${indels_recal%recal.vcf*}tranches

    local indels_tranches=${3:-$default_tranches}
    local indels_vqsr=${4:-$default_vqsr}

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


## add_snps_vqsr()
#
# Usage:
#	add_snps_vqsr sites_vcf snps_vcf snps_tranches snps_vqsr
#
# Apply SNP Recalibration
#
# Use the SNP recalibration model to add the SNP VQSR to the specified
# VCF file
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function add_snps_vqsr() {
    if [ $# -ne 4 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_vcf snps_vcf snps_tranches snps_vqsr"
        return
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3 $4"
    fi

    local sites_vcf=${1:-$joint_gvcf_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vqsr.vcf}
    local snps_recal=${2:-$joint_gvcf_dir/joint_gvfc.variant_filtered.sites_only.snps.recal.vcf}

    local default_vqsr=${sites_vcf%.vcf*}.snps.recal.vqsr.vcf
    local snps_tranches=${3:-$default_tranches}

    local default_tranches=${snps_recal/%recal.vcf/tranches}
    local snps_vqsr=${4:-$default_vqsr}

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


## collect_metrics()
#
# Usage:
#	collect_metrics sites_only_vcf intervals_list metrics_prefix
#
# collect variant calling metrics
#
#	Results will be stored in a set of files named starting with 
# 'metrics_prefix'
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function collect_metrics() {
    if [ $# -ne 3 ] ; then
        echo "Usage: ${FUNCNAME[0]} sites_only_vcf intervals_list metrics_prefix" ; return
    else
        echo ">>> ${FUNCNAME[0]} $1 $2 $3"
    fi

    local recalibrated_vcf=${1:-$joint_gvcf_dir/joint_gvfc.variant_filtered.sites_only.indels.recal.vqsr.snps.recal.vqsr.vcf}
    local intervals=$2

    local default_metrics=${1%vcf*}metrics
    local metrics_prefix=${3-$default_metrics}
    #metrics_prefix="$joint_gvcf_dir/joint_gvcfs_metrics"
    
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



## joint_genotyping()
#
# Usage:
#	joint_genotyping [gvcf_dir]
#
# Do a joint genotyping of all GVCF files present in the specified directory
#
#	If no directory is provided, a global default will be used instead,
# and if not global default is defined, then a directory named 'g.vcf' will
# be sought for VCF files.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
function joint_genotyping() {
    # process the command line
    if [ $# -gt 1 ] ; then
        echo "Usage: ${FUNCNAME[0]} [gvcf_dir]"
	return
    elif [ $# -eq 1 ] ; then
        gvcf_dir=$1
        echo ">>> ${FUNCNAME[0]} GVCF_DIR OVERRIDEN! Using $gvcf_dir"
        # do also override the names of the output directories
        joint_gvcf_dir="joint_${gvcf_dir}"
	joint_gvcf_out_dir="${joint_gvcf_dir}_analysis.mg=$MG"
    else
        gvcf_dir='g.vcf'
        echo ">>> ${FUNCNAME[0]} GVCF_DIR OVERRIDEN! Using $gvcf_dir"
        # do also override the names of the output directories
        joint_gvcf_dir="joint_${gvcf_dir}"
	joint_gvcf_out_dir="${joint_gvcf_dir}_analysis.mg=$MG"
    fi
    echo ">>> ${FUNCNAME[0]} $*"

    if [ ! -d "$gvcf_dir" ] ; then
    	echo ">>> ${FUNCNAME[0]} ERROR: '${gvcf_dir}' does not exist!"
        return
    fi

    if [ ! -d $joint_gvcf_dir ] ; then
        mkdir -p $joint_gvcf_dir
    fi
    if [ ! -d $joint_gvcf_out_dir ] ; then
        mkdir -p $joint_gvcf_out_dir
    fi
    
    workspace_dir_name=$joint_gvcf_dir

    # Make names by successively removing the trailing .vcf or .vcf.gz
    # and adding a new suffix. This will generate very long names, but
    # will also allow us to track the progress of the calculations.
    # ${var%vcf*} will remove the shortest match of 'vcf*' from the
    # end of the string (e.g. kk.vcf , kk.vcf.gz or kk.vcf-something 
    # will all yield kk.)

                                                      joint_gvcf=$joint_gvcf_out_dir/joint_gvcf.vcf.gz
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


    create_sample_name_map_from_GVCF_dir $gvcf_dir

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
    # NOTE: at this point some key annotation disappears (e.g. HaplotypeCaller added like GQ)
    make_sites_only_vcf $joint_gvcf_filtered_nocall \
		    $joint_gvcf_filtered_nocall_sites_only

    recalibrate_indels $joint_gvcf_filtered_sites_only \
		    $joint_gvcf_filtered_sites_only_indels

    recalibrate_snps $joint_gvcf_filtered_sites_only \
		    $joint_gvcf_filtered_sites_only_snps

## fails here with "no tranches above threshold"
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

}

#
#	FINALLY!
#
#	DO THE WORK
#
#	We do it this way for it allows us to work on the script while
# it is running, as this is the last line of the script.
#
joint_genotyping $*
