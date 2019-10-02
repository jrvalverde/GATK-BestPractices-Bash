#!/bin/bash
# 
# We use the GATK HaplotypeCaller to perform variant calling. The
# HaplotypeCaller is capable of calling SNPs and indels simultaneously via
# local de-novo assembly of haplotypes in an active region. In other words,
# whenever the program encounters a region showing signs of variation, it
# discards the existing mapping information and completely reassembles the
# reads in that region. This step is designed to maximize sensitivity in order
# to minimize false negatives, i.e. failing to identify real variants.

java -jar $GATK_JAR -T HaplotypeCaller \
	-R ${reference}.fna \
	-I $recal_reads.bam -o $raw_variants.vcf

# Extract SNPs and Indels
# 
# This step separates SNPs and Indels so they can be processed and analyzed
# independently.

# GATK3
java -jar $GATK_JAR -T SelectVariants \
       -R ${reference}.fna \
       -V $raw_variants.vcf \
       -selectType SNP \
       -o $raw_snps.vcf

java -jar $GATK_JAR -T SelectVariants \
	-R ${reference}.fna \
	-V $raw_variants.vcf \
	-selectType INDEL \
	-o $raw_indels.vcf

# GATK4
gatk -T SelectVariants \
	-R ${reference}.fna \
	-V $raw_variants.vcf \
	-selectType SNP \
	-o $raw_snps.vcf

gatk -T SelectVariants \
	-R ${reference}.fna \
	-V $raw_variants.vcf \
	-selectType INDEL \
	-o $raw_indels.vcf


# Filter Variants
# 
# The first step is designed to maximize sensitivity and is thus very lenient
# in calling variants. This is good because it minimizes the chance of missing
# real variants, but it means that we need to filter the raw callset produced
# in order to reduce the amount of false positives, which can be quite large.
# This is an important step in order to obtain the the highest-quality call set
# possible.
# 
# We apply the recommended hard filters for SNPs and Indels.

# GATK3
java -jar $GATK_JAR -T VariantFiltration \
	-R ${reference}.fna \
	-V $raw_snps.vcf \
	-filterName "QD_filter" \
	-filter "QD<2.0" \
	-filterName "FS_filter" \
	-filter "FS>60.0" \
	-filterName "MQ_filter" \
	-filter "MQ<40.0" \
	-filterName "SOR_filter" \
	-filter "SOR>10.0" \
	-o $filtered_snps.vcf

# Note: SNPs which are filtered out at this step will remain in the
# filtered_snps.vcf file, however, they will be marked as \u2018*_filter'
# based on which filter the SNP failed, while SNPs which passed the filter will
# be marked as PASS

java -jar $GATK_JAR -T VariantFiltration \
	-R ${reference}.fna \
	-V $raw_indels.vcf \
	-filterName "QD_filter" \
	-filter "QD<2.0" \
	-filterName "FS_filter" \
	-filter "FS>200.0" \
	-filterName "SOR_filter" \
	-filter "SOR>10.0" \
	-o $filtered_indels.vcf

# Note: Indels which are filtered out at this step will remain in
# the filtered_snps.vcf file, however, they will be marked as
# \u2018*_filter' based on which filter the indel failed, while Indels
# which passed the filter will be marked as PASS

# GATK4
gatk -T VariantFiltration \
	-R ${reference}.fna \
	-V $raw_snps.vcf \
	-filterName "QD_filter" \
	-filter "QD'<'2.0" \
	-filterName "FS_filter" \
	-filter "FS'>'60.0" \
	-filterName "MQ_filter" \
	-filter "MQ'<'40.0" \
	-filterName "SOR_filter" \
	-filter "SOR'>'10.0" \
	-o $filtered_snps.vcf

gatk -T VariantFiltration \
	-R ${reference}.fna \
	-V $raw_indels.vcf \
	-filterName "QD_filter" \
	-filter "QD'<'2.0" \
	-filterName "FS_filter" \
	-filter "FS'>'200.0" \
	-filterName "SOR_filter" \
	-filter "SOR'>'10.0" \
	-o $filtered_indels.vcf

# Annotation
# 
# We will use the program SnpEff. It annotates and predicts the effects of
# variants on genes (such as amino acid changes).
# 
# SnpEff has pre-built databases for thousands of genomes. We will be using the
# pre-built database GRCh38.p2.RefSeq.
# 
# Caveat:
# 
# Although we are using the NCBI RefSeq SnpEff Database, there is still an
# incompatibility in the chromosome names (ie. in our data, chromosome 20 is
# named NC_000020.11 whereas in the prebuilt SnpEff database, chromosome 20 is
# named '20'). To get around this, we will do a little magic to
# transform "NC_000020.11" to '20'.
# 
# First let's make a copy of filtered_snps.vcf:

# Locating and downloading the SnpEff database for your organism:
# 
# To locate and download the SnpEff database for your organism, execute the
# following command:

java -jar /share/apps/snpeff/4.1g/snpEff.jar databases \
    | grep -i 'Homo sapiens'

# Select the reference/version you are working with (for example:
# athalianaTair10), then:

java -jar snpEff.jar download $ref_snpeff_genome_refseq

cp $filtered_snps.vcf $filtered_snps_renamed.vcf

# We will edit our new file:
	
$EDIT $filtered_snps_renamed.vcf

# Use this snippet to replace all instanced of 'NC_000020.11' with '20' in vim:
# :%s/NC_000020.11/20/g
# :wq

# We're now ready to run SnpEff:
java -jar $snpEFF_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_snps_renamed.vcf \
	> $filtered_snps.ann.vcf

# Visualization (IGV)
# 
# Let's fire up IGV, load the Human_GRCh38 genome, and load our dedup.bam
# and filtered_snps.vcf files. Let's look at:
# 
# 1) Examples of high confidence SNPs in our filtered_snps.vcf file vs SNPs in
# our alignment not called by GATK
# 2) Examples of Homozygous variants vs Heterozygous variants.
# 3) Examples of SNPs filtered out by variant filtering (23, 703, 878)
# 4) Examples of SNPs in the cystatin C gene (CST3, variants in this gene are
# associated with Alzheimer's)

