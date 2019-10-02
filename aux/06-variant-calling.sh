set -ue
#
#       Variant Calling GATK 3.8
#       

proj='pr1'
#proj='pr2'
#proj='ref'

JAVAMEM=-Xms6000m
GATK3=~/contrib/gatk3
GATK4=~/contrib/gatk4
GATKVERS=3

KNDATA=~/work/lorena/data/gatk-hg19-bundle

genomeref=../data/ucsc.hg19.fasta
intervals=../panel/panel_65genes.ucsc.interval.list
outdir=variants-bwa
outdirra=$outdir/01-realign
outdirrc=$outdir/02-recalibration
#outdirug=$outdir/03-unifiedGenotyper
outdirhc=$outdir/03-haplotypeCaller
outdirvr=$outdir/04-variantRecalibrator
outdirhf=$outdir/05-HardFiltering

mkdir -p $outdirra
mkdir -p $outdirrc
#mkdir -p $outdirug
mkdir -p $outdirhc
mkdir -p $outdirvr
mkdir -p $outdirhf


for ali in bwa ; do
    for file in dedup-$ali/*ucsc*dedup.bam ; do
    
        outname=`basename $file .bam`

        if [ "$GATKVERS" -eq 3 ] ; then
            echo "Precomputing realignment for $file"
            if [ ! -e $outdirra/$outname.realign.list ] ; then
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T RealignerTargetCreator \
                    -R $genomeref \
                    -known $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -known $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
    	            -L $intervals \
                    -I $file \
                    -o $outdirra/$outname.realign.list
	        if [ $? -ne 0 ] ; then continue ; fi
	    fi

        else
           if [ ! -e $outdirra/$outname.realign.list ] ; then
               echo "Realigning is not needed in GATK4 if HaplotypeCaller or Mutect2 are used"
           
		echo $GATK4/gatk RealignerTargetCreator \
                    -R $genomeref \
                    -known $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -known $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
    	            -L $intervals \
                    -I $file \
                    -o $outdirra/$outname.realign.list IS UNSUPPORTED, USING VERSION 3 INSTEAD!!!
                    
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T RealignerTargetCreator \
                    -R $genomeref \
                    -known $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -known $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
    	            -L $intervals \
                    -I $file \
                    -o $outdirra/$outname.realign.list
	        if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi


        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirra/$outname.realign.bam ] ; then
                echo "Actually realigning BAM $file"
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T IndelRealigner \
                    -R $genomeref \
                    -I $file \
                    -known $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -known $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
                    -targetIntervals $outdirra/$outname.realign.list \
                    -o $outdirra/$outname.realign.bam
	        if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirrc/$outname.recal.grp ] ; then
                # pre-changes
                echo "Calibrating changes (precomputing)"
                java -jar $GATK3/GenomeAnalysisTK.jar \
                    -T BaseRecalibrator \
                    -R $genomeref \
                    -I $outdirra/$outname.realign.bam \
    	            -L $intervals \
                    -knownSites $KNDATA/dbsnp_138.hg19.vcf \
                    -knownSites $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -knownSites $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
                    -o $outdirrc/$outname.recal.grp 
                    
	        if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirrc/$outname.post_recal_data.grp ] ; then
                # post-changes
                echo "Calibrating changes (applying calibration)"
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T BaseRecalibrator \
                    -R $genomeref \
                    -I $outdirra/$outname.realign.bam \
    	            -L $intervals \
                    -knownSites $KNDATA/dbsnp_138.hg19.vcf \
                    -knownSites $KNDATA/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf \
                    -knownSites $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
                    -BQSR $outdirrc/$outname.recal.grp \
                    -o $outdirrc/$outname.post_recal_data.grp 
                    
	        if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirrc/$outname.realign.recal.bam ] ; then
                echo "Calibrating changes (printing new realigned recalibrated BAM)"
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T PrintReads \
                    -R $genomeref \
                    -I $outdirra/$outname.realign.bam \
    	            -L $intervals \
                    -BQSR $outdirrc/$outname.recal.grp \
                    -o $outdirrc/$outname.realign.recal.bam
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -lt 3 ] ; then
	    if [ ! -e $outdirrc/$outname.realign.recal.reduced.bam ] ; then
                echo "NOT DONE ON GATK >= 3"
	        java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T ReduceReads \
                    -R $genomeref \
                    -I $outdirrc/$outname.realign.recal.bam \
    	            -L $intervals \
                    -o $outdirrc/$outname.realign.recal.reduced.bam
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirhc/$outname.realign.raw_variants.vcf ] ; then
                echo "Calling variants on realigned file"
#                    -T HaplotypeCaller \
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                    -T UnifiedGenotyper \
                    -R $genomeref \
                    -I $outdirra/$outname.realign.bam \
                    -L $intervals \
                    -A Coverage \
                    -A QualByDepth \
                    --dbsnp $KNDATA/dbsnp_138.hg19.vcf \
                    -stand_call_conf 30 \
                    -o $outdirhc/$outname.realign.raw_variants.vcf
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirhc/$outname.realign.raw_variants.annotated.vcf ] ; then
                echo "Annotating variants of the realigned file"
                java -jar $GATK3/GenomeAnalysisTK.jar \
                  -T VariantAnnotator \
                  -R $genomeref \
                  -I $outdirra/$outname.realign.bam \
                  -V $outdirhc/$outname.realign.raw_variants.vcf \
                  -o $outdirhc/$outname.realign.raw_variants.annotated.vcf \
                  -A Coverage \
                  -A QualByDepth \
                  -all \
                  -L $outdirhc/$outname.realign.raw_variants.vcf \
                  --dbsnp $KNDATA/dbsnp_138.hg19.vcf
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirhc/$outname.realign.recal.raw_variants.vcf ] ; then
                echo "Calling variants of the realigned recalibrated file"
#                    -T HaplotypeCaller \
                java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
		    -T UnifiedGenotyper \
                    -R $genomeref \
                    -I $outdirrc/$outname.realign.recal.bam \
                    -L $intervals \
		    --dbsnp $KNDATA/dbsnp_138.hg19.vcf \
                    -stand_call_conf 30 \
                    -o $outdirhc/$outname.realign.recal.raw_variants.vcf \
                    -A Coverage \
                    -A QualByDepth 
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirhc/$outname.realign.recal.raw_variants.annotated.vcf ] ; then
                echo "Annotating variants of the realigned recalibrated file"
                java -jar $GATK3/GenomeAnalysisTK.jar \
                  -T VariantAnnotator \
                  -R $genomeref \
                  -I $outdirrc/$outname.realign.recal.bam \
                  -V $outdirhc/$outname.realign.recal.raw_variants.vcf \
                  -o $outdirhc/$outname.realign.recal.raw_variants.annotated.vcf \
                  -A Coverage \
                  -A QualByDepth \
                  -mvq 0.0 \
                  -L $outdirhc/$outname.realign.recal.raw_variants.vcf \
                  --dbsnp $KNDATA/dbsnp_138.hg19.vcf
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi

        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirvr/$outname.recalibrate_SNP.recal ] ; then
             java $JAVAMEM -jar $GATK3/GenomeAnalysisTK.jar \
                 -T VariantRecalibrator \
                 -R $genomeref \
                 -input $outdirhc/$outname.realign.recal.raw_variants.annotated.vcf \
                 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \
                     $KNDATA/hapmap_3.3.hg19.sites.vcf \
                 -resource:omni,known=false,training=true,truth=false,prior=12.0 \
                     $KNDATA/1000G_omni2.5.hg19.sites.vcf \
                 -resource:1000G,known=false,training=true,truth=false,prior=10.0 \
                     $KNDATA/1000G_phase1.indels.hg19.sites.vcf \
                 -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \
                     $KNDATA/dbsnp_138.hg19.vcf \
                 -an DP \
                 -an QD \
                 -an FS \
                 -an MQRankSum \
                 -an ReadPosRankSum \
                 -an QD -an MQ -an GQ \
                 -mode SNP \
                 -L $intervals \
		 -minNumBad 1000 \
                 -recalFile $outdirvr/$outname.recalibrate_SNP.recal \
                 -tranchesFile $outdirvr/$outname.recalibrate_SNP.tranches \
                 -rscriptFile $outdirvr/$outname.recalibrate_SNP_plots.R
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi


        if [ "$GATKVERS" -eq 3 ] ; then
	    if [ ! -e $outdirvr/$outname.recalibrated_snps_raw_indels.vcf ] ; then
            java -jar $GATK3/GenomeAnalysisTK.jar \
                -T ApplyRecalibration \
                -R $genomeref \
                -input $outdirhc/$outname.raw_variants.vcf \
                -mode SNP \
                --ts_filter_level 99.0 \
                -recalFile  $outdirvr/$outname.recalibrate_SNP.recal \
                -tranchesFile  $outdirvr/$outname.recalibrate_SNP.tranches \
                -o  $outdirvr/$outname.recalibrated_snps_raw_indels.vcf
	        
                if [ $? -ne 0 ] ; then continue ; fi
	    fi
        fi
exit

            java -jar GenomeAnalysisTK.jar \
                -T VariantRecalibrator \
                -R $genomeref \
                -input $outdirvr/$outname.recalibrated_snps_raw_indels.vcf \
                -resource:mills,known=true,training=true,truth=true,prior=12.0 mills.vcf \
                -an DP \
                -an FS \
                -an MQRankSum \
                -an ReadPosRankSum \
                -an QD -an MQ -an GQ \
                -mode INDEL \
                -tranche [100.0, 99.9, 99.0, 90.0] \
                -percentBad 0.01 \
                -minNumBad 1000 \
                -maxGaussians 4 \
                -recalFile $outdirvr/$outname.recalibrate_INDEL.recal \
                -tranchesFile $outdirvr/$outname.recalibrate_INDEL.tranches \
                -rscriptFile $outdirvr/$outname.recalibrate_INDEL_plots.R

            java -jar GenomeAnalysisTK.jar \
                -T ApplyRecalibration \
                -R $genomeref \
                -input $outdirvr/$outname.recalibrated_snps_raw_indels.vcf \
                -mode INDEL \
                --ts_filter_level 99.0 \
                -recalFile $outdirvr/$outname.recalibrate_INDEL.recal \
                -tranchesFile $outdirvr/$outname.recalibrate_INDEL.tranches \
                -o $outdirvr/$outname.recalibrated_variants.vcf

            java -jar GenomeAnalysisTK.jar \
                -T SelectVariants \
                -R $genomeref \
                -V $outdirhc/$outname.raw_variants.vcf \
                -L $intervals \
                -selectType SNP \
                -o $outdirhf/$outname.raw_snps.vcf

            java -jar GenomeAnalysisTK.jar \
                -T SelectVariants \
                -R $genomeref \
                -V $outdirhc/$outname.raw_HC_variants.vcf \
                -L $intervals \
                -selectType INDEL \
                -o $outdirhf/$outname.raw_indels.vcf

            java -jar GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $genomeref \
                -V $outdirhf/$outname.raw_snps.vcf \
                --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 ||  \
                HaplotypeScore > 13.0 || MappingQualityRankSum < -12.5 ||  \
                ReadPosRankSum < -8.0" \
                --filterName "my_snp_filter" \
                -o $outdirhf/$outname.filtered_snps.vcf

            java -jar GenomeAnalysisTK.jar \
                -T VariantFiltration \
                -R $genomeref \
                -V $outdirhf/$outname.raw_indels.vcf \
                --filterExpression "QD < 2.0 || FS > 200.0 || \
                ReadPosRankSum < -20.0" \
                --filterName "my_indel_filter" \
                -o $outdirhf/$outname.filtered_indels.vcf

    done
# Variables will take values from a configuration file. They will 
# depend on the experiment being analyzed
done

