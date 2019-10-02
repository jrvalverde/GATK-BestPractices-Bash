while read ch s e g ; do
    #echo $ch $s $e $g
    start=$s
    while [ $start -lt $e ] ; do
        end=$((start + 120))
        if [ $end -gt $e ] ; then end=$e ; fi
        echo -e "${ch}\t${start}\t${end}\t${g}"
        start=$end
    done
done
