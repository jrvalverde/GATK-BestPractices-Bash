mkdir tmp
GATK3=~/contrib/gatk3

for file in *dedup.bam ; do
    echo $file
    sample=`basename $file | sed -e 's/^*_S//g' -e 's/_.*//g'`
    # This should be used for alignment
    rginfo="@RG	ID:$proj	SM:$sample	PL:illumina	LB:$proj	PU:$proj"

    # Replace Read Group info with the one we want now, just in case
    # the alignment was not done with the info
    echo "Adding RG INFO: '$rginfo'"

    if [ ! -e tmp/$file ] ; then
        java -jar $GATK3/picard.jar AddOrReplaceReadGroups \
            VALIDATION_STRINGENCY=SILENT \
            I=$file \
            O=tmp/$file \
            SO=coordinate \
            RGID=1 \
            RGLB=proradium \
            RGPL=illumina \
            RGPU=unit1 \
            RGSM=$sample \
            CREATE_INDEX=TRUE
	if [ $? -ne 0 ] ; then echo "ERROR! HORROR!" ; exit ; fi
    fi
done

mv tmp/* .

