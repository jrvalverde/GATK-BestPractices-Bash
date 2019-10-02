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
 for i in vcf/*filtered.vcf ; do
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
#ref_snpeff_genome_RefSeq=hg19

sample=`echo $target | sed -e 's/.*_S/S/g' -e 's/_L*//g'`

vcf_dir=`dirname $target`.snpEff.ann.$ref_snpeff_genome_RefSeq
if [ ! -d $vcf_dir ] ; then mkdir $vcf_dir ; fi
cat=cat
# check if the file is compressed
ext=${target##*.}
if [ "$ext" == 'gz' ] ; then
    cat=zcat
    vcf_base=`basename $target .gz`
else
    vcf_base=`basename $target`
fi
# check that $vcf_base is a VCF file
ext=${vcf_base##*.}
if [ "$ext" != "vcf" ] ; then
    echo "ERROR: need a VCF file (.vcf or .vcf.gz) found $ext"
    exit
fi
    
vcf_base=$vcf_dir/`basename $vcf_base .vcf`	# remove extension

bwa=/usr/bin/bwa
PICARD_JAR=$HOME/contrib/gatk3/picard.jar
gatk=$HOME/contrib/gatk4/gatk
snpEff_path=$HOME/contrib/snpEff/snpEff

# Annotate VCF file
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
$cat $target | sed -e 's/chr//g' > $vcf_base.renamed.vcf

# We're now ready to run SnpEff:
if [ ! -s $vcf_base.ann.vcf ] ; then
    echo "java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$vcf_base.renamed.vcf \
	> $vcf_base.ann.vcf"
    java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$vcf_base.renamed.vcf \
	> $vcf_base.ann.vcf
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

