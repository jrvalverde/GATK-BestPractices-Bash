#!/bin/bash
# 
# We use the GATK HaplotypeCaller to perform variant calling. The
# HaplotypeCaller is capable of calling SNPs and indels simultaneously via
# local de-novo assembly of haplotypes in an active region. In other words,
# whenever the program encounters a region showing signs of variation, it
# discards the existing mapping information and completely reassembles the
# reads in that region. This step is designed to maximize sensitivity in order
# to minimize false negatives, i.e. failing to identify real variants.

target=$1
if [ "$target" == '' ] ; then
 for i in align/*.aligned.duplicates_marked.recalibrated.bam ; do
     echo $0 $i
     $0 $i
 done
 exit
fi


data=$HOME/work/lorena/data
gatk_bundle_dir=$data/gatk-hg19-bundle

reference=$data/ucsc.hg19
exome_intervals_list=exome.intervals.list
ref_snpeff_genome_RefSeq=GRCh37.p13.RefSeq
ref_snpeff_genome_RefSeq=hg19

sample=`echo $target | sed -e 's/.*_S/S/g' -e 's/_L*//g'`
read1=$target
read2=`echo $target | sed -e 's/_R1_/_R2_/g'`

align_base=align/`basename $target .bam`
aligned_reads=$align_base.aligned
sorted_reads=$aligned_reads.sorted
dedup_reads=$sorted_reads.dedup
metrics=$dedup_reads.metrics
recal_data=$dedup_reads.recal_data
recal_reads=$dedup_reads.recal

recal_reads=`dirname $target`/`basename $target .bam`
vcf_base='vcf'
raw_variants=$vcf_base/`basename $target .aligned.duplicates_marked.recalibrated.bam`.raw
raw_snps=$raw_variants.snps
raw_indels=$raw_variants.indels
filtered_snps=$raw_snps.filtered
filtered_indels=$raw_indels.filtered

bwa=/usr/bin/bwa
PICARD_JAR=$HOME/contrib/gatk3/picard.jar
gatk=$HOME/contrib/gatk4/gatk
snpEff_path=$HOME/contrib/snpEff/snpEff

if [ ! -d $vcf_base ] ; then mkdir -p vcf ; fi

if [ ! -s $raw_variants.vcf ] ; then
    $gatk HaplotypeCaller \
	-R ${reference}.fasta \
	-I $recal_reads.bam \
	-O $raw_variants.vcf \
    #	-ERC GVCF
fi
if [ ! -s $raw_variants.vcf ] ; then
    echo "Failed to make $raw_variants.vcf"
    exit
fi
# Extract SNPs and Indels
# 
# This step separates SNPs and Indels so they can be processed and analyzed
# independently.

if [ ! -s $raw_snps.vcf ] ; then
    $gatk SelectVariants \
	-R ${reference}.fasta \
	-V $raw_variants.vcf \
	--select-type SNP \
	-O $raw_snps.vcf
fi
if [ ! -s $raw_snps.vcf ] ; then 
    echo "Failed to make $raw_snps.vcf"
    exit
fi

if [ ! -s $raw_indels.vcf ] ; then
    $gatk SelectVariants \
	-R ${reference}.fasta \
	-V $raw_variants.vcf \
	-select-type INDEL \
	-O $raw_indels.vcf
fi
if [ ! -s $raw_indels.vcf ] ; then
    echo "Failed to make $raw_indels.vcf"
    exit
fi

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

# GATK4
if [ ! -s $filtered_snps.vcf ] ; then
    $gatk VariantFiltration \
	-R ${reference}.fasta \
	-V $raw_snps.vcf \
	--filter-name "QD_filter" \
	--filter-expression "QD < 2.0" \
	--filter-name "FS_filter" \
	--filter-expression "FS > 60.0" \
	--filter-name "MQ_filter" \
	--filter-expression "MQ < 40.0" \
	--filter-name "SOR_filter" \
	--filter-expression "SOR > 10.0" \
	-O $filtered_snps.vcf
fi
if [ ! -s $filtered_snps.vcf ] ; then
    echo "Failed to make $filtered_snps.vcf"
    exit
fi

# Note: SNPs which are filtered out at this step will remain in the
# filtered_snps.vcf file, however, they will be marked as \u2018*_filter'
# based on which filter the SNP failed, while SNPs which passed the filter will
# be marked as PASS

if [ ! -s $filtered_indels.vcf ] ; then
    $gatk VariantFiltration \
	-R ${reference}.fasta \
	-V $raw_indels.vcf \
	--filter-name "QD_filter" \
	--filter-expression "QD < 2.0" \
	--filter-name "FS_filter" \
	--filter-expression "FS > 200.0" \
	--filter-name "SOR_filter" \
	--filter-expression "SOR > 10.0" \
	-O $filtered_indels.vcf
fi
if [ ! -s $filtered_indels.vcf ] ; then
    echo "Failed to make $filtered_indels.vcf"
    exit
fi

# Note: Indels which are filtered out at this step will remain in
# the filtered_snps.vcf file, however, they will be marked as
# \u2018*_filter' based on which filter the indel failed, while Indels
# which passed the filter will be marked as PASS



# Annotation
# 
# We will use the program SnpEff. It annotates and predicts the effects of
# variants on genes (such as amino acid changes).
# 
# SnpEff has pre-built databases for thousands of genomes. We will be using the
# pre-built database GRCh38.p2.RefSeq.
# 

# Locating and downloading the SnpEff database for your organism:
# 
# To locate and download the SnpEff database for your organism, execute the
# following command:

#if [ ! -e snpEff_v4_3_$ref_snpeff_genome_RefSeq.zip ] ; then
if [ "YES" == "NO" ] ; then
    java -jar $snpEff_path/snpEff.jar databases \
        | grep -i 'Homo_sapiens'

# Select the reference/version you are working with (for example:
# athalianaTair10), GRCh37.p13.RefSeq, hg19, etc. then:

    echo java -jar $snpEff_path/snpEff.jar download $ref_snpeff_genome_RefSeq
    java -jar $snpEff_path/snpEff.jar download $ref_snpeff_genome_RefSeq

fi

# Caveat:
# 
# Although we are using the NCBI RefSeq SnpEff Database, there is still an
# incompatibility in the chromosome names (ie. in our data, chromosome 20 is
# named NC_000020.11 whereas in the prebuilt SnpEff database, chromosome 20 is
# named '20'). To get around this, we will do a little magic to
# transform "NC_000020.11" to '20'.
# 
# First let's make a copy of filtered_snps.vcf:
###cp $filtered_snps.vcf $filtered_snps.renamed.vcf

# We will edit our new file:
###$EDIT $filtered_snps.renamed.vcf

# Use this snippet to replace all instanced of 'NC_000020.11' with '20' in vim:
# :%s/NC_000020.11/20/g
# :wq

# in our case chromosomes are named 'chrNN' and should be NN
cat $filtered_snps.vcf | sed -e 's/chr//g' > $filtered_snps.renamed.vcf

# We're now ready to run SnpEff:
if [ ! -s $filtered_snps.ann.vcf ] ; then
    echo "java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_snps.renamed.vcf \
	> $filtered_snps.ann.vcf"
    java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_snps.renamed.vcf \
	> $filtered_snps.ann.vcf
fi

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

