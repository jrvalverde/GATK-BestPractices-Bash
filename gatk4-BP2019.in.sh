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
#
#	ABOUT THIS FILE:
#
# This file contains the configuration options for running the
# bash GATK best-practices from the Broad Institute scripts.
#
#	It is simply a BASH file. If you do not know much about BASH,
# then, please, read on.
#
#	It is line oriented.
#
#	On any line, everything after a '#' (hash-mark) is ignored and
# thus, can be used as a commentary.
#
#	A '\' symbol means 'take the next character directly' and not
# as its default meaning. Thus, if used before the end-of-a-line, that
# end-of-line will not be interpreted as such, hence effectively allowing
# you to continue it on the next line.
#
#	Values are assigned to variables using an '=' (equal sign).
#
#	Variables can be used by preceding them with a '$' (dollar sign)
# when there may be ambiguity, you should enclose the variable name in
# curly brackets ({ }).
#
#	When asigning a value that contains spaces, enclose it in double
# quotes (").
#
#	Some variables may take more than one value, these values are
# enclosed in parentheses '(' and ')' and separated between themselves
# by white space (e.g. ' ', TAB or newline).
#
#	If the same variable is assigned a value more than once, then the
# last value assigned will prevail
#
#

## SAMPLE NAME AND UNMAPPED BAMS
ref_name="Homo_sapiens_assembly38"
study_name="YOURSTUDYNAME"
prefix="$study_name"
library="YOURLIBRARYNAME"
sequencer="NAMEOFSEQUENCER"
date='2019-01-01T00:00:00+0100'
unmapped_bam_suffix=".bam"
flowcell_unmapped_bams_list="./uBam.list"

## PROGRAM-SPECIFIC PARAMETERS
SNP_VQSR_downsampleFactor=10
snp_filter_level=99.7
indel_filter_level=99.7

# We should only include here annotations that we want in the report
# `# comment` is a bash-ism to include comments in line.
# the DP annotation invoked by Coverage) should not be used when working 
# with exome datasets (gatk docs, args for VQSR)

#indel_recalibration_annotation_values=("FS" "ReadPosRankSum" "MQRankSum" "QD" "SOR" "DP")
# If any (e.g. MQRankSum) has zero variance, remove it
indel_recalibration_annotation_values=("FS" "ReadPosRankSum" "QD" "SOR" "DP")

#snp_recalibration_annotation_values=("QD" "MQRankSum" "ReadPosRankSum" "FS" "MQ" "SOR" "DP")
# If any (e.g. MQRankSum) has zero variance, remove it
snp_recalibration_annotation_values=("QD" "ReadPosRankSum" "FS" "MQ" "SOR" "DP")

recalibration_tranche_values=(
      100 99.95 99.9 99.8 99.7 99.6 99.5 99.4 99.3 99.0 98.0 97.0 05.0 90.0
      )


MG=2	# Max-gaussians parameter for VariantRecalibrator.
#MG=3   # The --max-gaussians parameter sets the expected number of clusters
#MG=4	# in modeling. If a dataset gives fewer distinct clusters, e.g. as can 
#MG=6   # happen for smaller data, then the tool will tell you there is 
        # insufficient data with a No data found error message. In this case, 
        # try decrementing the --max-gaussians value. 


## DATA DIRECTORIES
data_dir="./hg38"
gatk_bundle_dir="$data_dir/gatk_bundle"
broad_reference="$data_dir/broad-reference/v0"
panel_dir='./panel'
align_dir='align.gatk'
gvcf_dir='g.vcf'
joint_gvcf_dir="joint_${gvcf_dir}"
joint_gvcf_out_dir="${joint_gvcf_dir}_analysis.mg=$MG"


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


# default values for contamination data (use as reference NCBI Sequence)
#contamination_sites_ud="${VerifyBamID_home}/resource/hgdp.100k.b38.vcf.gz.dat.UD"
#contamination_sites_bed="${VerifyBamID_home}/resource/hgdp.100k.b38.vcf.gz.dat.bed"
#contamination_sites_mu="${VerifyBamID_home}/resource/hgdp.100k.b38.vcf.gz.dat.mu"
#contamination_sites_v="${VerifyBamID_home}/resource/hgdp.100k.b38.vcf.gz.dat.V"

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

omni_vcf="$gatk_bundle_dir/1000G_omni2.5.hg38.vcf.gz"

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
   
 
