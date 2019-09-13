set -ue
#
#       Check alignment and coverage metrics in ${ali} alignments
#	with and without adaptors
#
GATK=~/contrib/gatk3
GATK_VERS=3

#
#	FOR NOW WE PROCESS UCSC ALIGNMENTS WITH UCSC REFERENCE
#

for ali in bwa ; do
    OUT_DIR=dedup-${ali}-metrics
    if [ ! -d $OUT_DIR ] ; then mkdir $OUT_DIR ; fi

    for FILE in dedup-${ali}/*-ucsc-${ali}-srt-dedup.bam ; do
        BASE_NAME=`basename $FILE .bam`

        CASM_DIR=$OUT_DIR/CollectAlignmentSummaryMetrics
        mkdir -p $CASM_DIR
	CASM=$CASM_DIR/${BASE_NAME}_CollectAlignmentSummaryMetrics.txt
        if [ ! -e $CASM ] ; then
            echo ""
            echo "CollectAlignmentSummaryMetrics for $FILE"
            echo ""
 	    if [ "$GATK_VERS" -eq 4 ] ; then
                $GATK/gatk CollectAlignmentSummaryMetrics \
                    -R ../data/ucsc.hg19.fasta \
                    -I $FILE \
                    -O $CASM
            else
	        # picard.jar CollectAlignmentSummaryMetrics
                java -jar $GATK/picard.jar CollectAlignmentSummaryMetrics \
                    R=../data/ucsc.hg19.fasta \
                    I=$FILE \
                    O=$CASM
            fi
        fi

	CWM_DIR=$OUT_DIR/CollectWgsMetrics
        mkdir -p $CWM_DIR
        CWM=$CWM_DIR/${BASE_NAME}_CollectWgsMetrics.txt
        if [ ! -e $CWM ] ; then
            echo ""
            echo "Alignment Metrics with CollectWgsMetrics: $FILE  "
            echo ""
 	    if [ "$GATK_VERS" -eq 4 ] ; then
                gatk CollectWgsMetrics \
                    -I $FILE \
                    -O $CWM \
                    -R ../data/ucsc.hg19.fasta \
		    -INTERVALS ../panel/panel_65genes.ucsc.interval.list \
                    -INCLUDE_BQ_HISTOGRAM true
            else
                java -Xms6000m -jar $GATK/picard.jar \
                    CollectWgsMetrics \
                    I=$FILE \
                    O=$CWM  \
                    R=../data/ucsc.hg19.fasta \
                    INTERVALS= ../panel/panel_65genes.ucsc.interval.list \
                    INCLUDE_BQ_HISTOGRAM=true
            fi
        fi
        
        #echo Alignment Metrics with genomeCoverageBed: $FILE
        #genomeCoverageBed -ibam `basename $FILE` -g chromInfo_hg38_clean.txt > `basename $FILE_genomeCoverageBed.txt`
	SAM_DIR=$OUT_DIR/samdepth
        mkdir -p $SAM_DIR
        SAM=$SAM_DIR/${BASE_NAME}.samdepth
        if [ ! -e $SAM ] ; then
	    echo ""
            echo Alignment Metrics with samtools: $FILE
            echo ""
                     
            samtools depth $FILE \
                > $SAM
                
            cat $SAM \
                | awk '$3 >= 20 {print ;}' \
                > $SAM_DIR/$BASE_NAME.samdepth.20+
            cat $SAM \
                | awk '$3 >= 100 {print ;}' \
                > $SAM_DIR/$BASE_NAME.samdepth.100+

            #cat $SAM \
            #    | awk '{sum+=$3} END { print "Tumor Average = ",sum/NR}' \
            #    > $SAM_DIR/$BASE_NAME.samdepth.coverage


        fi

        SMB_DIR=$OUT_DIR/sambamba
        mkdir -p $SMB_DIR
        SBB=$SMB_DIR/${BASE_NAME}_sambamba.txt
        continue
        if [ ! -e $SBB ] ; then
            echo ""
            echo Alignment Metrics with sambamba: $FILE 
            echo ""

            sambamba depth base \
                    -L ../panel/panel_65genes.bed \
        	    -F 'mapping_quality >= 20' -t 4 -c 0 \
                    -o $SBB $FILE \
		    #-m	# fix mate overlaps
        fi

    done

done
