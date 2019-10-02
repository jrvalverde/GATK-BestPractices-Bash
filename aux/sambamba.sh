FILE=$1
BASE_NAME=`basename $FILE .bam`

PATH=$HOME/work/lorena/bin:$PATH

echo "------------------------------------------------------------------"
echo " Computing coverage depths of $FILE"
echo "------------------------------------------------------------------"


OUT_DIR=stats
if [ ! -d $OUT_DIR ] ; then mkdir -p $OUT_DIR ; fi


SAM_DIR=$OUT_DIR/samdepth
mkdir -p $SAM_DIR
if [ ! -s $SAM_DIR/$BASE_NAME.samdepth.100+ ] ; then
    #echo Alignment Metrics with genomeCoverageBed: $FILE
    #genomeCoverageBed -ibam `basename $FILE` -g chromInfo_hg38_clean.txt > `basename $FILE_genomeCoverageBed.txt`
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

        cat $SAM \
            | awk '{sum+=$3} END { print "Tumor Average = ",sum/NR}' \
            > $SAM_DIR/$BASE_NAME.samdepth.coverage


    fi
fi


SMB_DIR=$OUT_DIR/sambamba
mkdir -p $SMB_DIR
SBB=$SMB_DIR/${BASE_NAME}_sambamba.txt
if [ ! -e $SBB ] ; then
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
fi
