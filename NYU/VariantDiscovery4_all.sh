#!/bin/bash
# 
# We use the GATK HaplotypeCaller to perform variant calling. The
# HaplotypeCaller is capable of calling SNPs and indels simultaneously via
# local de-novo assembly of haplotypes in an active region. In other words,
# whenever the program encounters a region showing signs of variation, it
# discards the existing mapping information and completely reassembles the
# reads in that region. This step is designed to maximize sensitivity in order
# to minimize false negatives, i.e. failing to identify real variants.

if [ $# -eq 0 ] ; then
    $0 align/*.recal*.bam
    exit
elif [ $# -gt 1 ] ; then
    for i in  "$@" ; do
        $0 $i
    done
    exit
fi
# $# must be exactly 1

target=$1

data="./hg19"
gatk_bundle_dir="$data/gatk_bundle"

reference="$data/ucsc.hg19"
ref_dict="$data/ucsc.hg19.dict"
exome_intervals_list="exome.intervals.list"
exome_picard_list="exome.ucsc.hg19.interval_list"
ref_snpeff_genome_RefSeq="hg19"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"

sample=`echo $target | sed -e 's/.*_S/S/g' -e 's/_L.*//g'`

recal_reads=`dirname $target`/`basename $target .bam`
vcf_base='vcf'
raw_variants=$vcf_base/`basename $target .aligned.sorted.dedup.recal.bam`.raw
raw_snps=$raw_variants.snps
raw_indels=$raw_variants.indels
filtered_snps=$raw_snps.filtered
filtered_indels=$raw_indels.filtered
annotated_snps=$filtered_snps.ann
annotated_indels=$filtered_indels.ann

gatk=$HOME/contrib/gatk4/gatk
snpEff_path=$HOME/contrib/snpEff/snpEff

if [ ! -d $vcf_base ] ; then mkdir -p vcf ; fi

if [ ! -s $raw_variants.vcf ] ; then
    $gatk HaplotypeCaller \
	-R ${reference}.fasta \
	-I $recal_reads.bam \
	-O $raw_variants.vcf \
        --dbsnp $dbSNP_vcf \
        --emit-ref-confidence GVCF \
        #-L $exome_picard_list	# if commented we get all
fi
if [ ! -s $raw_variants.vcf ] ; then
    echo "Failed to make $raw_variants.vcf"
    exit
fi


############################################################################
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


############################################################################
# Filter Variants
# 
# The first step is designed to maximize sensitivity and is thus very lenient
# in calling variants. This is good because it minimizes the chance of missing
# real variants, but it means that we need to filter the raw callset produced
# in order to reduce the amount of false positives, which can be quite large.
# This is an important step in order to obtain the highest-quality call set
# possible.
# 
# We first apply the recommended hard filters for SNPs and Indels.

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
# filtered_snps.vcf file, however, they will be marked as '*_filter'
# based on which filter the SNP failed, while SNPs which passed the 
# filter will be marked as PASS

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
# '*_filter' based on which filter the indel failed, while Indels
# which passed the filter will be marked as PASS


# Final filtration step to add putative LPZ filters
### LPZ
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
if [ ! -s "$filtered_indels.lpz.vcf" ] ; then
    $gatk VariantFiltration \
	-R ${reference}.fasta \
	-V "$filtered_indels.vcf" \
        --filter-name "HighHaplotypeScore" \
        --filter-expression "HaplotypeScore > 13.0" \
        --filter-name "HighQualZero" \
        --filter-expression "(MQ0 / DP) > 0.1" \
        --filter-name "IndelFilter" \
        --filter-expression "QD < 2.0 || FS > 200.0 || ReadPosRankSum < -20.0" \
        --filter-name "LowGQ" \
        --filter-expression "GQ < 20" \
        --filter-name "LowMQRankSum" \
        --filter-expression "MQRankSum < -12.5" \
        --filter-name "LowQualByDepth" \
        --filter-expression 'QD < 2.0' \
        --filter-name "LowRMSMappingQual" \
        --filter-expression 'MQ < 30.0' \
        --filter-name "LowReadDepth" \
        --filter-expression 'DP < 16' \
        --filter-name "LowReadPosRankSum" \
        --filter-expression 'ReadPosRankSum < -8.0' \
        --filter-name "StrandBias" \
        --filter-expression 'FS > 60.0' \
	-O "$filtered_indels.lpz.vcf"
fi
if [ ! -s "$filtered_indels.lpz.vcf" ] ; then
    echo "Failed to make $filtered_indels.lpz.vcf"
    exit
fi



############################################################################
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
cat $filtered_indels.vcf | sed -e 's/chr//g' > $filtered_indels.renamed.vcf
cat $filtered_indels.lpz.vcf | sed -e 's/chr//g' > $filtered_indels.lpz.renamed.vcf

# We're now ready to run SnpEff on SNPs:
if [ ! -s $annotated_snps.vcf ] ; then
    echo "java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_snps.renamed.vcf \
	> $annotated_snps.vcf"
    java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_snps.renamed.vcf \
	> $annotated_snps.vcf
fi
if [ ! -s $filtered_snps.ann.vcf ] ; then
    echo "Failed to make $filtered_snps.ann.vcf"
    exit
fi

# Now annotate indels
if [ ! -s $annotated_indels.vcf ] ; then
    echo "java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_indels.renamed.vcf \
	> $annotated_indels.vcf"
    java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_indels.renamed.vcf \
	> $annotated_indels.vcf
fi
if [ ! -s $filtered_indels.ann.vcf ] ; then
    echo "Failed to make $filtered_indels.ann.vcf"
    exit
fi

# And finally annotate ultra-filtered indels:
if [ ! -s $filtered_indels.lpz.ann.vcf ] ; then
    echo "java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_indels.lpz.renamed.vcf \
	> $filtered_indels.lpz.ann.vcf"
    java -jar $snpEff_path/snpEff.jar \
	-v $ref_snpeff_genome_RefSeq \
	$filtered_indels.lpz.renamed.vcf \
	> $filtered_indels.lpz.ann.vcf
fi
if [ ! -s $filtered_indels.lpz.ann.vcf ] ; then
    echo "Failed to make $filtered_indels.lpz.ann.vcf"
    exit
fi


############################################################################
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



############################################################################
# Variant metric calculation
#
for input_vcf in ${raw_variants}*.ann.vcf ; do
    base=`echo $input_vcf | sed -e 's/\.vcf$//g'`
    # make index if missing (annotate does not make one)
    if [ ! -e $input_vcf.idx ] ; then
        $gatk IndexFeatureFile -F $input_vcf
    fi
    if [ ! -s ${base}.variant_calling_summary_metrics \
      -o ! -s ${base}.variant_calling_detail_metrics ] ; then
	# for GATK v4
        $gatk --java-options "-Xms2000m" \
            CollectVariantCallingMetrics \
            --INPUT $input_vcf \
            --OUTPUT $base \
            --DBSNP $dbSNP_vcf \
            --TARGET_INTERVALS ${exome_picard_list} \
            --REFERENCE_SEQUENCE $reference.fasta \
            --GVCF_INPUT true \
            #--SEQUENCE_DICTIONARY $ref_dict 
    fi
    if [ ! -s ${base}.variant_calling_summary_metrics \
      -o ! -s ${base}.variant_calling_detail_metrics ] ; then
        echo "Couldn't collect metrics for $base"
        exit
    fi
done
