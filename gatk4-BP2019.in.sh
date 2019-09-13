## SAMPLE NAME AND UNMAPPED BAMS
ref_name="ucsc.hg19"
study_name="PR2"
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

#MG=6
#MG=4
#MG=3
MG=2	# Max-gaussians parameter for VariantRecalibrator.
        # The --max-gaussians parameter sets the expected number of clusters
	# in modeling. If a dataset gives fewer distinct clusters, e.g. as can 
        # happen for smaller data, then the tool will tell you there is 
        # insufficient data with a No data found error message. In this case, 
        # try decrementing the --max-gaussians value. 


## DATA DIRECTORIES
data_dir="./hg19"
gatk_bundle_dir="$data_dir/gatk_bundle"
panel_dir='./panel'
align_dir='align.gatk'
gvcf_dir='g.vcf'
workspace_dir_name='joint_gvcf'
out_dir="joint_gvcf_analysis.mg=$MG"


## PATHS
bwa_path='/usr/bin/'
gatk_exec="${HOME}/contrib/gatk4/gatk"
gotc_path="${HOME}/contrib/gatk4/"
gitc_path="${HOME}/contrib/gatk3/"


## AUXILIARY FILES
# may be a Picard interval_list, a GATK interval or list, or a BED file
# but a Picard interval_list is recommended
study_interval_list="./panel/panel5-120.ucsc.hg19.interval_list"
#
# This one needs to be in GATK intervals/list format so 
# we can add unmapped at the end
study_intervals="./panel/panel5-120.ucsc.hg19.intervals"


## "KNOWN SITES RESOURCES"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"
dbSNP_vcf_index="$gatk_bundle_dir/dbsnp_138.hg19.vcf.idx"
hapmap_vcf="$gatk_bundle_dir/hapmap_3.3.hg19.sites.vcf"
omni_vcf="$gatk_bundle_dir/1000G_omni2.5.hg19.sites.vcf"
Mills_vcf="$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
G1000_vcf="$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"

known_indels_sites_VCFs=(
    "$gatk_bundle_dir/dbsnp_138.hg19.vcf"
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"
#	Homo_sapiens_assembly38.known_indels.vcf.gz	N/A for hg19
)
known_indels_sites_indices=(
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf.idx"
#	Homo_sapiens_assembly38.known_indels.vcf.gz.tbi	N/A for hg19
)
   
