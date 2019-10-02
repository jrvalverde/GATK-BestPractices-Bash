while read ch s e g ; do
    #echo $ch $s $e $g
    start=$s
    while [ $start -lt $e ] ; do
        end=$((start + 120))
        if [ $end -ge $e ] ; then end=$e ; fi
        # For bed files we need zero-offset coordinates
        coord1=$((start - 1))
        coord2=$((end - 2))
        #echo -e "${ch}\t${start}\t${end}\t${g}"
        echo -e "${ch}\t${coord1}\t${coord2}\t${g}"
        start=$end
    done
done
