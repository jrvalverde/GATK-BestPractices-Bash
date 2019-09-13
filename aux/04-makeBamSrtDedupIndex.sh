set -ue

GATK=~/contrib/gatk

GATK_VERS=3

align_dir=align

if [ ! -d $align_dir ] ; then mkdir -p $align_dir ; fi


for ali in bwa ; do
    if [ ! -d $align_dir ] ; then 
        echo "No $ali alignments to process" 
        continue
    fi
    if [ ! -d dedup-$ali ] ; then mkdir dedup-$ali ; fi
    for FILE in $align_dir/*.sam ; do
        echo "Sorting/Deduplicating/Indexing $FILE"
        SRT=dedup-$ali/`basename $FILE .sam`-srt.bam
        DEDUP=dedup-$ali/`basename $SRT .bam`-dedup.bam
        BAI=dedup-$ali/`basename $DEDUP .bam`.bai

        if [ $GATK_VERS -eq 3 ] ; then
	    if [ ! -e "$SRT" ] ; then
                #Sort
                java -Xms6g -jar $GATK/picard.jar SortSam \
                	INPUT=$FILE \
                        OUTPUT="$SRT" \
                        SORT_ORDER=coordinate
	    fi

            if [ ! -e "$DEDUP" ] ; then
                #Mark and remove duplicates
                java -Xms6g -jar $GATK/picard.jar MarkDuplicates \
                	REMOVE_DUPLICATES=true \
                        INPUT="$SRT" \
                        OUTPUT="$DEDUP" \
                        METRICS_FILE=metrics.txt
            fi

            if [ ! -e "$BAI" ] ; then
                #Indexing
                java -Xms6g -jar $GATK/picard.jar BuildBamIndex \
                	INPUT="$DEDUP"
	    fi

        elif [ $GATK_VERS -eq 4 ] ; then   
	    if [ ! -e "$SRT" ] ; then
                #Sort
                echo "Sorting"
                $GATK/gatk SortSam \
                	-I $FILE \
                        -O $SRT \
                        -SORT_ORDER coordinate
            fi

            if [ ! -e "$DEDUP" ] ; then
                #Mark and remove duplicates
                echo "Deduplicating"
                $GATK/gatk MarkDuplicates \
                	-REMOVE_DUPLICATES true \
                        -I "$SRT" \
                        -O "$DEDUP" \
                        --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            		-METRICS_FILE metrics.txt
            fi

            if [ ! -e "$BAI" ] ; then
                #Indexing
                echo "Indexing"
                $GATK/gatk BuildBamIndex -I "$DEDUP"
            fi
        fi
    done
done
