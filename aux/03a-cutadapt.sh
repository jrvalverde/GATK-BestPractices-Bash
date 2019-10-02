#
#	Example script to remove adaptors by using cutadapt sofware
#	This step needs to be performed before the alignment
#
#	Could we include tipical sequence adaptors used from illumina 
#	in case we are not sure about our adaptor sequence?
#	Only used to test if this step improves the alignment metrics or if 
#	it could affect the variant calling
#


set -ue
minoverlap=3
minlength=20
for FILE in 18106FL-01-02-14_S31_L002 18106FL-01-02-15_S32_L002 18106FL-01-02-16_S33_L002 18106FL-01-02-17_S34_L002 18106FL-01-03-01_S3_L003 18106FL-01-03-02_S36_L003 18106FL-01-03-03_S37_L003 18106FL-01-03-04_S38_L003 18106FL-01-03-05_S39_L003 18106FL-01-03-06_S40_L003 
do
echo Remove adapters from $FILE
cutadapt -j 0 -m $minlength --match-read-wildcards -O $minoverlap -b AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -B AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o "$FILE"_R1_001_wo_adapt_l20_ok.fastq.gz -p "$FILE"_R2_001_wo_adapt_l20_ok.fastq.gz "$FILE"_R1_001.fastq.gz "$FILE"_R2_001.fastq.gz > "$FILE"_cutadapt.log
done
