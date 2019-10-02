#!/bin/bash

bamfile=$1
if [ "$bamfile" == "" ] ; then
    for i in align/*.recalibrated.bam ; do
        $0 $i
    done
    exit
fi
set -ue

seqdb=$HOME/work/lorena/data/seqdb.hg19
RefSeq_locdb=$HOME/work/lorena/data/RefSeq.locdb
ref_genome_fasta=$HOME/work/lorena/data/ucsc.hg19.fasta

gatk3_dir=$HOME/contrib/gatk3
xhmm_dir=cnv.xhmm

export PATH=~/work/lorena/bin:$PATH

# to ease our work
if [ ! -d $xhmm_dir ] ; then
    mkdir $xhmm_dir
fi
if [ ! -e EXOME.interval.list ] ; then
    ln -s ../panel/INGEM_panel5_120_mod_interval_format.intervals EXOME.interval.list
fi

out_base=$xhmm_dir/`basename $bamfile .bam`

# Run GATK to calculate depth of coverage for each sample
echo ">>> calculating depth coverage for $bamfile"
if [ ! -s $out_base.DATA.sample_statistics ] ; then
    echo $bamfile
    echo $out_base

    java -Xmx4g -jar $gatk3_dir/GenomeAnalysisTK.jar \
            -T DepthOfCoverage \
            -I "$bamfile" \
            -L EXOME.interval.list \
            -R $ref_genome_fasta \
            -dt BY_SAMPLE \
            -dcov 300 \
            -l INFO \
            --omitDepthOutputAtEachBase \
            --omitLocusTable \
            --minBaseQuality 0 \
            --minMappingQuality 20 \
            --start 1 --stop 5000 \
            --nBins 200 \
            --includeRefNSites \
            --countType COUNT_FRAGMENTS \
            -o "$out_base".DATA \
#           -mmq 30 -mbq 30 -dels \
#           -ct 1 -ct 10 -ct 20 -ct 30 -ct 40 -ct 50 \
#           -omitLocusTable
fi

# Now that we have the output in a separate directory we can work inside it
cd $xhmm_dir

if [ ! -e ./EXOME.interval.list ] ; then
    ln -s ../EXOME.interval.list .
fi

# we'll use the reference seqdb from the Harvard PLINK/SEQ site in
# https://atgu.mgh.harvard.edu/plinkseq/resources.shtml
# http://psychgen.u.hpc.mssm.edu/plinkseq_resources/hg19/seqdb.hg19.gz
if [ ! -e seqdb ] ; then
    ln -s $seqdb seqdb
fi

# Same for locdb
# http://psychgen.u.hpc.mssm.edu/plinkseq_resources/hg19/locdb.gz
if [ ! -e RefSeq.locdb ] ; then
    ln -s $RefSeq_locdb RefSeq.locdb
fi

if [ ! -e ucsc.hg19.fasta ] ; then
    ln -s $ref_genome_fasta ucsc.hg19.fasta
fi

if [ ! -e ~/work/lorena/scripts/example_make_XHMM_plots.R ] ; then
    ln -s ~/work/lorena/scripts/example_make_XHMM_plots.R .
fi

# parameters for CNV
# 
# XHMM parameters file
# A parameters file consists of the following 9 values:
# 
#     Exome-wide CNV rate
#     Mean number of targets in CNV
#     Mean distance between targets within CNV (in KB)
#     Mean of DELETION z-score distribution
#     Standard deviation of DELETION z-score distribution
#     Mean of DIPLOID z-score distribution
#     Standard deviation of DIPLOID z-score distribution
#     Mean of DUPLICATION z-score distribution
#     Standard deviation of DUPLICATION z-score distribution
# 
# As an example, the file with parameters:
# 
# ***********************************************************************
# Input CNV parameters file:
# ***********************************************************************
# 1e-08   6       70      -3      1       0       1       3       1
# ***********************************************************************
# 
# translates into XHMM parameters of:
# 
# ***********************************************************************
# Pr(start DEL) = Pr(start DUP) = 1e-08
# Mean number of targets in CNV [geometric distribution] = 6
# Mean distance between targets within CNV [exponential decay] = 70 KB
# 
# DEL read depth distribution ~ N(mean=-3, var=1)
# DIP read depth distribution ~ N(mean=0, var=1)
# DUP read depth distribution ~ N(mean=3, var=1)
# ***********************************************************************
# 
cat > params.txt <<END
1e-8	6	70	-3	1.00	0	1.00	3	1.00
END

### jr ###
cat > params.txt <<END
1e-8	200	70	-3	1.00	0	1.00	3	1.00
END


# Use xhmm to combine GATK Depth-of-coverage outputs for multiple samples
# (at the same loci)

echo ">>> merging samples"
if [ ! -s DATA.RD.txt ] ; then
    # prepare list of files to merge
    #args=''
    echo -n ''  > DATA.sample_interval_summary.files.txt
    for i in *.DATA.sample_interval_summary ; do
        #args="$args -GATKdepths $i"
        #args="$arg $i"
        echo "$i" >> DATA.sample_interval_summary.files.txt
    done

    #xhmm --mergeGATKdepths -o xhmm-out/DATA.RD.txt $args
    xhmm  -o DATA.RD.txt --mergeGATKdepths \
        --GATKdepthsList=DATA.sample_interval_summary.files.txt
fi

# Run GATK to calculate the per-target GC contents and create a list of
# the targets with extreme GC content

echo ">>> analysing GC content"
if [ ! -e extreme_gc_targets.txt ] ; then
    java -Xmx3072m -jar $gatk3_dir/GenomeAnalysisTK.jar \
        -T GCContentByInterval \
        -L EXOME.interval.list \
        -R $ref_genome_fasta \
        -o DATA.locus_GC.txt

    cat DATA.locus_GC.txt \
        | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' \
        > extreme_gc_targets.txt
fi

# Run PLink/Seq to calculate the fraction of repeat-masked bases in each
# target and create a list of those to filter out

echo ">>> analysing repeat-masked bases"
if [ ! -e EXOME.targets.reg ] ; then
    ~/work/lorena/scripts/interval_list_to_pseq_reg \
        EXOME.interval.list > EXOME.targets.reg
fi

if [ ! -s ./EXOME.targets.LOCDB.loc-load ] ; then
    pseq . loc-load --locdb ./EXOME.targets.LOCDB \
        --file EXOME.targets.reg --group targets \
        --out ./EXOME.targets.LOCDB.loc-load

fi

if [ ! -s low.complexity.targets.txt ] ; then
    pseq . loc-stats --locdb ./EXOME.targets.LOCDB \
        --group targets --seqdb ./seqdb \
        | awk '{if (NR > 1) print $_ }' \
        | sort -k1 -g \
        | awk '{print $10}' \
        | paste ./EXOME.interval.list - \
        | awk '{print $1"\t"$2}' \
        > DATA.locus_complexity.txt

    cat DATA.locus_complexity.txt \
        | awk '{if ($2 > 0.25) print $1}' \
        > low.complexity.targets.txt
fi

# Filter samples and targets and then mean-center the targets

echo ">>> filtering samples and targets"
if [ -s low.complexity.targets.txt ] ; then
  if [ ! -s ./DATA.filtered_centered.RD.txt ] ; then
    xhmm --matrix -r DATA.RD.txt \
        --centerData --centerType target \
        -o ./DATA.filtered_centered.RD.txt \
        --outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
        --excludeTargets ./extreme_gc_targets.txt \
        --excludeTargets ./low.complexity.targets.txt \
        --minTargetSize 10 --maxTargetSize 10000 \
        --minMeanTargetRD 10 --maxMeanTargetRD 500 \
        --minMeanSampleRD 25 --maxMeanSampleRD 200 \
        --maxSdSampleRD 150
  fi
  # In our case, filtering produces no data, so we'll avoid
  # filtering, plus I don't know the correct parameters to use
    xhmm --matrix -r DATA.RD.txt \
        --centerData --centerType target \
        -o ./DATA.filtered_centered.RD.txt \
        --outputExcludedTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
        --outputExcludedSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \

fi

# Run PCA on mean-centered data

echo ">>> running PCA"
if [ ! -s DATA.RD_PCA ] ; then
    xhmm --PCA -r DATA.filtered_centered.RD.txt --PCAfiles ./DATA.RD_PCA

fi

# Normalize mean-centered data using PCA information

echo ">>> normalizing"
if [ ! -e DATA.PCA_normalized.txt ] ; then
    xhmm --normalize \
        -r ./DATA.filtered_centered.RD.txt \
        --PCAfiles ./DATA.RD_PCA \
        --normalizeOutput DATA.PCA_normalized.txt \
        --PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7
fi

# Filter and Z-score center (by sample) the PCA-normalized data

echo ">>> filtring and Z-score centering"
if [ ! -s DATA.PCA_normalized.filtered.sample_zscores.RD.txt ] ; then
    xhmm --matrix \
        -r DATA.PCA_normalized.txt \
        --centerData --centerType sample --xScoreData \
        -o DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
        --outputExcludedTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
	--outputExcludedSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
	--maxSdTargetRD 30
fi

# Filter original read-depth data to be the same as filtered, normalized data

echo ">>> filtering original data"
if [ ! -s DATA.same_filtered.RD.txt ] ; then
    xhmm --matrix -r ./DATA.RD.txt \
	--excludeTargets ./DATA.filtered_centered.RD.txt.filtered_targets.txt \
	--excludeTargets ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
	--excludeSamples ./DATA.filtered_centered.RD.txt.filtered_samples.txt \
	--excludeSamples ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
	-o ./DATA.same_filtered.RD.txt
fi

# Discover CNVs in normalized data

echo ">>> discovering CNVs"
if [ ! -s DATA.xcnv ] ; then
    xhmm --discover -p ./params.txt \
	-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
	-R ./DATA.same_filtered.RD.txt \
	-c ./DATA.xcnv \
	-a ./DATA.aux_xcnv \
	-s ./DATA
fi

#XHMM .xcnv format - explained
#-----------------------------
# Consider the output from the tutorial example data:
# 
# SAMPLE    CNV  INTERVAL               KB      CHR   MID_BP     TARGETS     NUM_TARG   Q_EXACT   Q_SOME   Q_NON_DIPLOID   Q_START   Q_STOP    MEAN_RD   MEAN_ORIG_RD
# HG00121   DEL  22:18898402-18913235   14.83   22    18905818   104..117    14         9         90       90              8         4         -2.51     37.99
# HG00113   DUP  22:17071768-17073440   1.67    22    17072604   4..11       8          25        99       99              53        25        4.00      197.73
# 
# SAMPLE 	sample name
# CNV 		type of copy number variation (DEL or DUP)
# INTERVAL 	genomic range of the called CNV
# KB 		length in kilobases of called CNV
# CHR 		chromosome name on which CNV falls
# MID_BP 	the midpoint of the CNV (to have one genomic number for 
#               plotting a single point, if desired)
# TARGETS 	the range of the target indices over which the CNV is called 
#               (NOTE: considering only the FINAL set of post-filtering targets)
# NUM_TARG 	# of exome targets of the CNV
# Q_EXACT 	Phred-scaled quality of the exact CNV event along the entire 
#               interval
#               - Identical to EQ in .vcf output from genotyping
# Q_SOME 	Phred-scaled quality of some CNV event in the interval
#               - Identical to SQ in .vcf output from genotyping
# Q_NON_DIPLOID 	Phred-scaled quality of not being diploid, i.e., DEL 
#               or DUP event in the interval
#               - Identical to NDQ in .vcf output from genotyping
# Q_START 	Phred-scaled quality of "left" breakpoint of CNV
#               - Identical to LQ in .vcf output from genotyping
# Q_STOP 	Phred-scaled quality of "right" breakpoint of CNV
#               - Identical to RQ in .vcf output from genotyping
# MEAN_RD 	Mean normalized read depth (z-score) over interval
#               - Identical to RD in .vcf output from genotyping
# MEAN_ORIG_RD 	Mean read depth (# of reads) over interval
#               - Identical to ORD in .vcf output from genotyping


# Genotype discovered CNVs in all samples:

echo ">>> genotyping CNVs"
if [ ! -s DATA.vcf ] ; then
    xhmm --genotype -p ./params.txt \
	-r ./DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R ./DATA.same_filtered.RD.txt \
	-g ./DATA.xcnv -F ./ucsc.hg19.fasta \
	-v ./DATA.vcf
fi



# And the R visualization plots
#
# Annotate exome targets with their corresponding genes:

echo ">>> anotating targets"
if [ ! -e ./annotated_targets.refseq ] ; then
    pseq . loc-intersect \
	--group refseq \
	--locdb ./RefSeq.locdb \
	--file ./EXOME.interval.list \
	--out ./annotated_targets.refseq

fi

# Plot the XHMM pipeline and CNV discovered:

Rscript example_make_XHMM_plots.R



# Visualizing XHMM results with included R scripts
# ------------------------------------------------
# 
# The sources/scripts/ directory contains R scripts for visualizing the PCA
# results and the read depth data at each called CNV.
# 
# See the example_make_XHMM_plots.R script for an easy starting point to
# running these scripts.
# 
# The first plot you'll get is: mean_sample_coverage.pdf
# 
# This image gives the sample-wide distribution of exome-wide coverage, where
# each per-sample coverage value is the mean of the coverage values calculated
# for each exome target (which itself is the mean coverage at all of its bases
# in that particular sample). In this experiment, we sequenced each sample to a
# mean coverage of 150x, so that we expect a typical sample to indeed have 150
# reads covering an average base in an average exome target.
# 
# And, you'll also similarly get: mean_target_coverage.pdf
# 
# This analogously gives the target-wide distribution of coverage (over all
# samples). That is, each per-target coverage value is the mean of the
# per-sample coverage values at that target (where again, this is the mean
# coverage at all of its bases in that sample). As above, since our goal was to
# have 150x coverage exome-wide, we'd expect each target to have around 150x
# coverage, but we see here that there is high variability in target coverage.
# For example, some targets have as much 400x coverage (averaged over all
# samples), and we also see a non-trivial number of targets that have 0
# coverage for all samples (e.g., targets where capture has presumably failed).
# 
# Next, we consider the principal component analysis (PCA) normalization stage.
# We compare each of the principal components to known sample and target
# features: PC_correlations.pdf
# 
# The dotted line (at PC = 15) indicates that XHMM automatically removed the
# 1st 15 components based on their significant relative variance. Now, in this
# plot we consider known sample and target features (that XHMM did not
# incorporate in its decision to remove them). We see that these 1st 15 PC tend
# to show correlation with various target features (colored circles) such as GC
# content and the mean depth of sequencing coverage at that target, and also
# with various sample features (colored diamonds) such as gender and mean depth
# of sequencing for that sample. On the other hand, there is a marked change in
# quality of the PC after the 1st 15 or so, with a sudden drop-off in the
# levels of correlation with genome-wide and batch effects expected to strongly
# bias the read depth of coverage.
# 
# Related to this is the "scree" plot for the PCA, showing the standard
# deviation of the depth data independently ascribed to each of the principal
# components: PC_stddev.pdf
# 
# This case is typical, where we see that the cut-off automatically detected by
# XHMM corresponds to a significant drop in the variance (an "elbow" in the
# curve). Note the log scale of the y axis.
# 
# In the output directory, the "PC" sub-directory is created, which contains
# plots of the read depth data projected into each of the principal components:
# PC/PC.*.png
# 
# This principal component (the 3rd one, in this case) has found the variance
# in read depth due to gender differences, with males having lower coverage on
# the X chromosome and higher coverage on the Y. Therefore, the loadings for
# this component have a correlation of 0.99 with the gender of the samples.
# 
# Next, before z-scores are calculated and the HMM is run to call CNV for each
# sample, we perform a final filtering step. Here, we remove any targets that
# have "very scattered" read depth distributions post normalization. These can
# be thought of as targets for which the normalization may have failed, and it
# is better to remove such strong effects (still likely to be artifacts) to
# prevent them from drowning out other more subtle signals. To do this, we
# remove any targets with large standard deviations of their post-normalization
# read depths across all samples. As a (proto-typical) example, we see here
# that the small fraction of targets with standard deviations any larger than
# the 30 to 50 range should be removed (in this case, for a scenario of ~ 100x
# mean sequencing coverage). This plot will be output to: per_target_sd.pdf
# 
# 
# 
# Lastly, in the "plot_CNV" sub-directory, a plot is made focusing on each CNV
# called in each individual (as found in the .xcnv file): plot_CNV/sample_*.png
# 
# The example above shows a de novo deletion called by XHMM (in red) spanning
# the DLGAP1 gene (discussed in this study and previously validated in this
# study).
# 
# In general, these plots show each sample's read depths at each of the targets
# in focus, which are connected by gray lines. If the sample has a called
# deletion, then it's colored in red, and duplications in green. Gene names are
# added below to annotate the genomic region.
# 
# Also, it's important to note that, by default, a plot is generated for EACH
# CNV called by XHMM, without considering its quality score or frequency in the
# sample. As for many analyses, it is incumbent upon the researcher to consider
# the metrics and properties of the data set as a whole before considering any
# one particular CNV call. That is, here we recommend filtering on both the SQ
# quality threshold (the "Q_SOME" column in the .xcnv file) and filtering on
# sample-wide frequency to obtain a high-quality call set. See the publications
# for details.
# 




# Finally return to previous folder
cd -



