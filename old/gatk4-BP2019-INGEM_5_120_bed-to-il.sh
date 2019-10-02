GATK4=$HOME/contrib/gatk4/

$GATK4/gatk BedToIntervalList \
	-I INGEM_panel5_120_mod.bed \
        -O INGEM_panel5_120_mod.ucsc.hg19.interval_list \
        -SD ../data/ucsc.hg19.fasta \
	-R ../data/ucsc.hg19.fasta
