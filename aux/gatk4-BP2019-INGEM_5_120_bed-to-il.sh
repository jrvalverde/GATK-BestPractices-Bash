GATK4=$HOME/contrib/gatk4/

bedfile=$1
ref_fasta=$2
output=${1%.bed}.interval_list

$GATK4/gatk BedToIntervalList \
	-I $bedfile \
        -O $output \
        -SD $ref_fasta \
	-R $ref_fasta
