## Copyright CNB-CSIC, 2019
## Copyright Broad Institute, 2019
## 
## This pipeline implements data pre-processing according to the GATK 
## Best Practices.  
##
## Requirements/expectations :
## - Pair-end sequencing data in unmapped BAM (uBAM) format
## - One or more read groups, one per uBAM file, all belonging to a single 
##	sample (SM)
## - Input uBAM files must additionally comply with the following 
## 	requirements:
##   	- filenames all have the same suffix (we use ".unmapped.bam")
##   	- files must pass validation by ValidateSamFile 
##   	- reads are provided in query-sorted order
##   	- all reads must have an RG tag
##
## Output :
## - A clean BAM file and its index, suitable for variant discovery analyses.
##
##
## LICENSING : 
## This script is released under the WDL source code license (BSD-3) (see LICENSE in 
## https://github.com/broadinstitute/wdl). Note however that the programs it calls may 
## be subject to different licenses. Users are responsible for checking that they are
## authorized to run all programs before running this script. Please see the dockers
## for detailed licensing information pertaining to the included programs.

DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

## SAMPLE NAME AND UNMAPPED BAMS
ref_name=human.genome.19
sample_name="patient.data"
flowcell_unmapped_bams_list="~/work/exome/patient.data/uBam.list"
unmapped_bam_suffix=".bam"
data_dir="/home/scientific/work/exome/data/"
gatk_bundle_dir="$data_dir/gatk-hg19-bundle"
bwa_path='/usr/bin/'
gatk_exec="/home/scientific/contrib/gatk4/gatk"
gotc_path="/home/scientific/contrib/gatk4/"


# get definitions from local file
if [ -s processing-for-variant-discovery-gatk4.in.sh ] ; then
    . processing-for-variant-discovery-gatk4.in.sh
else
    cp `dirname $0`/gatk4-BP2019-processing-for-variant-discovery.in.sh \
       ./processing-for-variant-discovery-gatk4.in.sh
    echo "I have created a preferences file in this directory:"
    echo ""
    echo "        processing-for-variant-discovery-gatk4.in.sh"
    echo ""
    echo "Please, open it in a text editor, adapt it to your needs and"
    echo "run this command again."
    exit
fi

## REFERENCE FILES
ref_dict="$data_dir/ucsc.hg19.dict"
ref_fasta="$data_dir/ucsc.hg19.fasta"
ref_fasta_index="$data_dir/ucsc.hg19.fasta.fai"
ref_sa="$data_dir/ucsc.hg19.fasta.sa"
ref_amb="$data_dir/ucsc.hg19.fasta.amb"
ref_bwt="$data_dir/ucsc.hg19.fasta.bwt"
ref_ann="$data_dir/ucsc.hg19.fasta.ann"
ref_pac="$data_dir/ucsc.hg19.fasta.pac"

## "KNOWN SITES RESOURCES"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"
dbSNP_vcf_index="$gatk_bundle_dir/dbsnp_138.hg19.vcf.idx"

known_indels_sites_VCFs=(
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"
#	Homo_sapiens_assembly38.known_indels.vcf.gz	N/A for hg19
)
known_indels_sites_indices=(
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.idx"
    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf.idx"
#	Homo_sapiens_assembly38.known_indels.vcf.gz.tbi	N/A for hg19
)
   
## "MISC PROGRAM PARAMETERS", 
compression_level=5
num_cpu="16"

## "PATHS",
align_bwa_dir="./align_bwa"


## "JAVA OPTIONS", 
SamToFastqAndBwaMem_java_opt="-Xms3000m"
MergeBamAlignment_java_opt="-Xms3000m"
MarkDuplicates_java_opt="-Xms4000m"
SortAndFixTags_java_opt_sort="-Xms4000m"
SortAndFixTags_java_opt_fix="-Xms500m"
BaseRecalibrator_java_opt="-Xms4000m"
GatherBqsrReports_java_opt="-Xms3000m"
ApplyBQSR_java_opt="-Xms3000m"
GatherBamFiles_java_opt="-Xms2000m"
GatherBamFiles_java_opt="-Xms8000m"
# default
java_opt="-Xms4000m"

## "MEMORY ALLOCATION", 
GetBwaVersion_mem_size="1 GB"
SamToFastqAndBwaMem_mem_size="14 GB"
MergeBamAlignment_mem_size="3500 MB"
MarkDuplicates_mem_size="7 GB"
SortAndFixTags_mem_size="5000 MB"
CreateSequenceGroupingTSV_mem_size="2 GB"
BaseRecalibrator_mem_size="6 GB"
GatherBqsrReports_mem_size="3500 MB"
ApplyBQSR_mem_size="3500 MB"
GatherBamFiles_mem_size="3 GB"
# default
mem_size="7 GB"

## ALL TOGETHER FILE NAME
base_file_name=$sample_name.$ref_name

## LIST OF BAM FILES TO PROCESS
#flowcell_unmapped_bams=$(< $flowcell_unmapped_bams_list)
# bash 4+
readarray -t flowcell_unmapped_bams < $flowcell_unmapped_bams_list
# note that we may as well read the file line by line instead of using an array


## TASK DEFINITIONS

# Get version of BWA
function GetBwaVersion () {
  
    # Not setting "set -o pipefail" here because /bwa has a rc=1 and we don't want to allow rc=1 to succeed 
    # because the sed may also fail with that error and that is something we actually want to fail on.
    ${bwa_path}/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
}


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment
function SamToFastqAndBwaMem () {

    set -o pipefail
    set -e
    
    # input_bam is a full-length absolute path name
    unmapped_bam=$1
    
    # Get the basename, i.e. strip the filepath and the extension
    bam_basename=`basename $unmapped_bam $unmapped_bam_suffix`
    output_unmerged_bam=$align_bwa_dir/$bam_basename.unmerged.bam
    
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

    bwa_commandline="bwa mem -p -v 3 -t ${num_cpu} -Y ${bash_ref_fasta}"

    if [ ! -s ${output_unmerged_bam} ] ; then
	mkdir -p $align_bwa_dir
        if [ ! -s $bam_basename.fastq ] ; then
	    echo "Generating $bam_basename.fastq"
            
            # convert to FastQ, pipe to bwa for alignment and to sam for bam comversion
            java -Dsamjdk.compression_level=${compression_level} \
    	          ${SamToFastqAndBwaMem_java_opt} \
                  -jar ${gotc_path}/picard.jar \
                  SamToFastq \
	          INPUT=${unmapped_bam} \
	          FASTQ=$bam_basename.fastq \
	          INTERLEAVE=true \
	          NON_PF=true 
	fi
        echo "Aligning and making ${output_unmerged_bam}"
        ${bwa_path}${bwa_commandline} $bam_basename.fastq   \
            2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
            > $bam_basename.sam
        samtools view -1 $bam_basename.sam > ${output_unmerged_bam}
        rm -f $bam_basename.fastq $bam_basename.sam
    fi

}


# Merge original input uBAM file with BWA-aligned BAM file
function MergeBamAlignment () {
    unmapped_bam=$1

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    bam_basename=`basename $unmapped_bam $unmapped_bam_suffix`
    unmerged_aligned_bam=$align_bwa_dir/$bam_basename.unmerged.bam
    output_bam=$align_bwa_dir/${bam_basename}.bam
    
    if [ ! -s $output_bam ] ; then
	echo "Merging BAM alignment with original unmapped bam"
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${MergeBamAlignment_java_opt}" \
          MergeBamAlignment \
          --VALIDATION_STRINGENCY SILENT \
          --EXPECTED_ORIENTATIONS FR \
          --ATTRIBUTES_TO_RETAIN X0 \
          --ALIGNED_BAM ${unmerged_aligned_bam} \
          --UNMAPPED_BAM ${unmapped_bam} \
          --OUTPUT ${output_bam} \
          --REFERENCE_SEQUENCE ${ref_fasta} \
          --PAIRED_RUN true \
          --SORT_ORDER "unsorted" \
          --IS_BISULFITE_SEQUENCE false \
          --ALIGNED_READS_ONLY false \
          --CLIP_ADAPTERS false \
          --MAX_RECORDS_IN_RAM 2000000 \
          --ADD_MATE_CIGAR true \
          --MAX_INSERTIONS_OR_DELETIONS -1 \
          --PRIMARY_ALIGNMENT_STRATEGY MostDistant \
          --PROGRAM_RECORD_ID "bwamem" \
          --PROGRAM_GROUP_VERSION "${bwa_version}" \
          --PROGRAM_GROUP_COMMAND_LINE "${bwa_commandline}" \
          --PROGRAM_GROUP_NAME "bwamem" \
          --UNMAPPED_READ_STRATEGY COPY_TO_TAG \
          --ALIGNER_PROPER_PAIR_FLAGS true \
          --UNMAP_CONTAMINANT_READS true
    fi
    # we no longer need the unmerged file
    ###rm $unmerged_aligned_bam
}

# Mark duplicate reads to avoid counting non-independent observations
# May work with only one or an array of many files
function MarkDuplicates() {
    
    input_bams=( $@ )
    
    if [ ${#input_bams[@]} -gt 1 ] ; then
        # use a summary name
        #base_file_name=$sample_name.$ref_name	# defined at the top
        output_bam_basename=${base_file_name}
    else
        # use a file specific name
        output_bam_basename=`basename $input_bams .bam`
    fi
    output_bam=$align_bwa_dir/$output_bam_basename.aligned.unsorted.duplicates_marked.bam
    metrics_filename=$align_bwa_dir/${output_bam_basename}.duplicate_metrics

    # prepare argument list of input files
    input_list=''
    for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done
    
    
    if [ ! -s ${output_bam} ] ; then
        # Task is assuming query-sorted input so that the Secondary and 
        # Supplementary reads get marked correctly. 
        # This works because the output of BWA is query-grouped and therefore,
        # so is the output of MergeBamAlignment. 
        # While query-grouped isn't actually query-sorted, it's good enough for
        # MarkDuplicates with ASSUME_SORT_ORDER="queryname"
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${MarkDuplicates_java_opt}" \
            MarkDuplicates \
            ${input_list} \
            --OUTPUT ${output_bam} \
            --METRICS_FILE ${metrics_filename} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true
    fi    
}


function MarkDuplicatesOld () {
    
    input_bams=( $@ )
    
    if [ ${#input_bams[@]} -gt 1 ] ; then
        # use a summary name
        #base_file_name=$sample_name.$ref_name	# defined at the top
        output_bam_basename=${base_file_name}
        output_bam=$align_bwa_dir/$output_bam_basename.aligned.unsorted.duplicates_marked.bam
        metrics_filename=$align_bwa_dir/${output_bam_basename}.duplicate_metrics

        # prepare argument list of input files
        input_list=''
        for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done
    else
        # use a file specific name
        output_bam_basename=`basename $input_bams .bam`
        output_bam=$align_bwa_dir/$output_bam_basename.aligned.unsorted.duplicates_marked.bam
        metrics_filename=$align_bwa_dir/${output_bam_basename}.duplicate_metrics

	input_list=' --INPUT $input_bams'
    fi
    
    
    if [ ! -s ${output_bam} ] ; then
        # Task is assuming query-sorted input so that the Secondary and 
        # Supplementary reads get marked correctly. 
        # This works because the output of BWA is query-grouped and therefore,
        # so is the output of MergeBamAlignment. 
        # While query-grouped isn't actually query-sorted, it's good enough for
        # MarkDuplicates with ASSUME_SORT_ORDER="queryname"
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${MarkDuplicates_java_opt}" \
            MarkDuplicates \
            ${input_list} \
            --OUTPUT ${output_bam} \
            --METRICS_FILE ${metrics_filename} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true
    fi
    
    return
    
    # That was for all bam files together (where each read has been properly
    # tagged. This is for each bam file separately (where each bam corresponds
    # to a different sample)
    #
    # This SHOULD NOT BE DONE HERE, but outside by repeatedly calling us
    #
    if [ ${#input_bams[@]} -gt 1 ] ; then
      for i in ${input_bams[@]} ; do 
        out=$align_bwa_dir/`basename $i .bam`.aligned.unsorted.duplicates_marked.bam
        metrics=align_bwa_dir/`basename $i .bam`.duplicate_metrics
        echo "Marking duplicates in $i"
        echo "Saving results in $out"
        if [ ! -s $out ] ; then
            ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${MarkDuplicates_java_opt}" \
                MarkDuplicates \
                --INPUT ${i} \
                --OUTPUT ${out} \
                --METRICS_FILE ${metrics} \
                --VALIDATION_STRINGENCY SILENT \
                --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
                --ASSUME_SORT_ORDER "queryname" \
                --CREATE_MD5_FILE true
        fi
      done
    fi
}


# Sort BAM file by coordinate order and fix tag values for NM and UQ
function SortAndFixTags() {

    input_bam=$1
    
    unsorted_bam_base=`basename $input_bam .aligned.unsorted.duplicates_marked.bam`
    output_bam=$align_bwa_dir/${unsorted_bam_base}.aligned.duplicate_marked.sorted.bam
  
    set -o pipefail

    #echo "Sorting >>> $input_bam <<<"
    if [ ! -s $output_bam ] ; then
        ${gatk_exec} \
          --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_sort}" \
          SortSam \
          --INPUT ${input_bam} \
          --OUTPUT /dev/stdout \
          --SORT_ORDER "coordinate" \
          --CREATE_INDEX false \
          --CREATE_MD5_FILE false \
        | \
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_fix}" \
          SetNmMdAndUqTags \
          --INPUT /dev/stdin \
          --OUTPUT ${output_bam} \
          --CREATE_INDEX true \
          --CREATE_MD5_FILE true \
          --REFERENCE_SEQUENCE ${ref_fasta}
    fi
    echo "Created $output_bam"  
}


# Generate sets of intervals for scatter-gathering over chromosomes
function CreateSequenceGroupingTSV() {

  if [ ! -s sequence_grouping.txt -o ! -s sequence_grouping_with_unmapped.txt ] 
  then
    # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
    # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
    # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    /home/scientific/contrib/anaconda2/bin/python <<CODE

with open("${ref_dict}", "r") as ref_dict_file:
    sequence_tuple_list = []
    longest_sequence = 0
    for line in ref_dict_file:
        if line.startswith("@SQ"):
            line_split = line.split("\t")
            # (Sequence_Name, Sequence_Length)
            sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
    longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
# We are adding this to the intervals because hg38 has contigs named with embedded colons (:) and a bug in 
# some versions of GATK strips off the last element after a colon, so we add this as a sacrificial element.
hg19_protection_tag = ":1+"
# initialize the tsv string with the first sequence
tsv_string = sequence_tuple_list[0][0] + hg19_protection_tag
temp_size = sequence_tuple_list[0][1]
for sequence_tuple in sequence_tuple_list[1:]:
    if temp_size + sequence_tuple[1] <= longest_sequence:
        temp_size += sequence_tuple[1]
        tsv_string += "\t" + sequence_tuple[0] + hg19_protection_tag
    else:
        tsv_string += "\n" + sequence_tuple[0] + hg19_protection_tag
        temp_size = sequence_tuple[1]
# add the unmapped sequences as a separate line to ensure that they are recalibrated as well
with open("sequence_grouping.txt","w") as tsv_file:
  tsv_file.write(tsv_string)
  tsv_file.close()

tsv_string += '\n' + "unmapped"

with open("sequence_grouping_with_unmapped.txt","w") as tsv_file_with_unmapped:
  tsv_file_with_unmapped.write(tsv_string)
  tsv_file_with_unmapped.close()

CODE
  fi

}

# Generate Base Quality Score Recalibration (BQSR) model
function BaseRecalibrator() {
    input_bam=$1
    sequence_group_interval=$2

    # either we do it for all sequence group intervals at once
    # or we generate one different recalibration_report for each
    # interval, so we can later gather all together in 
    # a call to GatherBsqrReports()
    #base_out_name=$base_file_name
    base_out_name=`basename $input_bam .aligned.duplicate_marked.sorted.bam`
    recalibration_report_filename=$align_bwa_dir/${base_out_name}.recal_data.csv
    known_indel_sites_arg=''
    for i in ${known_indels_sites_VCFs[@]} ; do
        known_indel_sites_arg="${known_indel_sites_arg} --known-sites $i"
    done
    sg_arg=''
    for i in ${sequence_group_interval[@]} ; do
        sg_arg="${sg_arg} -L $i"
    done
    
    if [ ! -s ${recalibration_report_filename} ] ; then
        ${gatk_exec} --java-options "${BaseRecalibrator_java_opt}" \
          BaseRecalibrator \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${recalibration_report_filename} \
          --known-sites ${dbSNP_vcf} \
          $known_indel_sites_arg \
          -L sequence_grouping.intervals.list
          #$sg_arg
    fi
}


# Combine multiple recalibration tables from scattered BaseRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked differently.
function GatherBqsrReports {
    input_bqsr_reports=( $@ )
  

    inrep_arg=''
    for i in ${input_bqsr_reports[@]} ; do
        inrep_arg="${inrep_arg} -I $i"
    done
    output_report_filename=$align_bwa_dir/${base_file_name}.all.recal_data.csv

    if [ ! -s ${output_report_filename} ] ; then
        ${gatk_exec} --java-options "${GatherBqsrReports_java_opt}" \
            GatherBQSRReports \
            $inrep_arg \
            -O ${output_report_filename}
    fi
}


# Apply Base Quality Score Recalibration (BQSR) model
# using a previously computed BQSR report (by BaseRecalibrator /
# GatherBqsrReports)

function ApplyBQSR() {
    input_bam=$1
    sequence_group_interval=$2

    # either we do it for all sequence group intervals at once
    # or we generate one different recalibration_report for each
    # interval, so we can later gather all together in 
    # a call to GatherBsqrReports()
    #base_out_name=$base_file_name
    base_name=`basename $input_bam .aligned.duplicate_marked.sorted.bam`
    recalibration_report=$align_bwa_dir/$base_name.recal_data.csv
    output_bam=$align_bwa_dir/$base_name.aligned.duplicates_marked.recalibrated.bam

    sg_arg=''
    for i in ${sequence_group_interval[@]} ; do
        sg_arg="${sg_arg} -L $i"
    done
    
    if [ ! -s ${output_bam} ] ; then
        ${gatk_exec} --java-options "${ApplyBQSR_java_opt}" \
             ApplyBQSR \
             -R ${ref_fasta} \
             -I ${input_bam} \
             -O ${output_bam} \
             $sg_arg \
             -bqsr ${recalibration_report} \
             --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 \
             --add-output-sam-program-record \
             --create-output-bam-md5 \
             --use-original-qualities
    fi
}


# Combine multiple recalibrated BAM files from scattered ApplyRecalibration 
# runs
function GatherBamFiles() {
    input_bams=( $@ )
  
    in_bams_arg=''
    for i in ${input_bams[@]} ; do  
        in_bams_arg="${in_bams_arg} --INPUT $i"
    done
    output_bam=./$align_bwa/$base_file_name.all.align.dup.recal.bam

    if [ ! -s $output_bam ] ; then
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${GatherBamFiles_java_opt}" \
             GatherBamFiles \
             ${in_bams_arg} \
             --OUTPUT ${output_bam} \
             --CREATE_INDEX true \
             --CREATE_MD5_FILE true
    fi
}


# WORKFLOW DEFINITION 

function PreProcessingForVariantDiscovery_GATK4() {
    # Get the version of BWA to include in the PG record in the header of the BAM produced 
    # by MergeBamAlignment. 
    bwa_version=`GetBwaVersion`

    # Align flowcell-level unmapped input bams
    #for unmapped_bam in "${flowcell_unmapped_bams[@]}" ; do
    for unmapped_bam in $(cat $flowcell_unmapped_bams_list) ; do

	echo "Processing $unmapped_bam"

        # Map reads to reference
        $CALL SamToFastqAndBwaMem $unmapped_bam

        # Merge original uBAM and BWA-aligned BAM 
        $CALL MergeBamAlignment $unmapped_bam
        
    done
    # get list of bam files generated into an array
    # we have generated a file "X.unmerged.bam" and an "X.bam" in the align dir
    merged_bams=$( ls $align_bwa_dir/*001.bam ) 
    
    # mark duplicates on a per-file basis
    for i in ${merged_bams[@]} ; do
        out=$align_bwa_dir/`basename $i .bam`.aligned.unsorted.duplicates_marked.bam
	if [ ! -s $out ] ; then
            $CALL MarkDuplicates $i
	fi
    done
    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    # We take advantage of the tool's ability to take multiple BAM inputs and 
    # write out a single output to avoid having to spend time just 
    # merging BAM files.
    $CALL MarkDuplicates "${merged_bams[@]}"
    marked_all_bam=$align_bwa_dir/${base_file_name}.aligned.unsorted.duplicates_marked.bam 



    # Sort deduped BAM files and fix tags
    # This will also consider the aggregate file
    for i in $align_bwa_dir/*.aligned.unsorted.duplicates_marked.bam ; do
        unsorted_bam_base=`basename $i .aligned.unsorted.duplicates_marked.bam`
        sorted_bam=$align_bwa_dir/$unsorted_bam_base.aligned.duplicate_marked.sorted.bam
        if [ ! -e $sorted_bam ] ; then
            $CALL SortAndFixTags $i
        fi
    done
    # the sorted-fixed aggregate bam file
    sorted_all_bam=$align_bwa_dir/${base_file_name}.aligned.duplicate_marked.sorted.bam

    # ------------------ ATTENTION -----------------
    # Create list of sequences for scatter-gather parallelization 
    CreateSequenceGroupingTSV 
    # read groupings into arrays
    sequence_grouping=$(< "sequence_grouping.txt")
    sequence_grouping_with_unmapped=$(< "sequence_grouping_with_unmapped.txt")

    # Here we may want to use our own groupings only for the genes we are
    # interested in, saving the global groupings and generating new ones:
    # We will group by chromosome
    USE_EXOME_COORDINATES='yes'
    if [ "$USE_EXOME_COORDINATES" == "yes" ] ; then
        if [ ! -e sequence_grouping.txt.sav ] ; then
            mv sequence_grouping.txt sequence_grouping.txt.sav
            mv sequence_grouping_with_unmapped.txt sequence_grouping_with_unmapped.txt.sav
            echo -n "" > sequence_grouping.txt
            cut -d: -f1 EXOME.intervals.list | sort | uniq |\
                while read chr ; do
                    echo "grouping $chr"
                    grep $chr EXOME.intervals.list | \
                    while IFS=$'\n' read interval ; do
                        echo -n "$interval	" >> sequence_grouping.txt
                    done
                    echo "">> sequence_grouping.txt
                done
            cp sequence_grouping.txt sequence_grouping_with_unmapped.txt
            echo "unmapped" >> sequence_grouping_with_unmapped.txt
        fi
    fi
    ### NOTE THAT WE'LL IGNORE THIS PART
    BY_SEQ_GRPS='NO'
    
    # ------------------ ATTENTION -----------------    
    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM
    # Perhaps we should do them all at once? Check documentation
    if [ "$BY_SEQ_GRPS" == 'YES' ] ; then
        echo "BaseRecalibrator by groups NOT YET READY"
        exit
        # THIS DOES NOT WORK: THE TOOL REFUSES TO TAKE THE SUBGROUPS 
        # GENERATED BY GATK (for chromosomes it should be chrX, not
        # chrX:1+ and for other regions, it is supposed to be
        # chrX:nnn-NNN). Yet it complains.
        # PLUS WE SHOULD SAVE EACH REPORT WITH A SEPARATE NAME
        #for subgroup in ${sequence_grouping[@]} ; do
        cat sequence_grouping.txt | while read subgroup ; do
          # Generate the recalibration model by interval
          $CALL BaseRecalibrator $sorted_all_bam "'$subgroup'"

          for sorted_bam in $align_bwa_dir/*.duplicate_marked.sorted.bam ; do
              $CALL BaseRecalibrator $sorted_bam "'$subgroup'"
          done
        done
    else
        echo "Using EXOME.intervals.list"
        if [ ! -s sequence_grouping.intervals.list ] ; then
            cp EXOME.intervals.list sequence_grouping.intervals.list
        fi
        # Generate the recalibration model by interval
        $CALL BaseRecalibrator $sorted_all_bam sequence_grouping.intervals.list
        for sorted_bam in $align_bwa_dir/*.duplicate_marked.sorted.bam ; do
            $CALL BaseRecalibrator $sorted_bam sequence_grouping.intervals.list
        done
    fi
    # get report name for full report (all samples together)
    base_report_name=`basename $sorted_all_bam .aligned.duplicate_marked.sorted.bam`
    recalibration_report=$align_bwa_dir/${base_report_name}.recal_data.csv

    # Merge the individual recalibration reports resulting from 
    # recalibration.
    # If we recalibrated by sequence groups, each one will produce a
    # recalibration report. If we did by individuals, recalibrating
    # each with all SGs at once, then each individual will produce a
    # recalibration reports.
    # We have
    #	- a recalibration report for all individuals merged with all SG
    #   - a report for each individual with all SG at once
    #	We'll gather the later individual reports into one reort
    individual_reports=$( ls $align_bwa_dir/*.recal_data.csv | grep -v $base_file_name )
    
    $CALL GatherBqsrReports "${individual_reports[@]}"
    output_bqsr_report=$align_bwa_dir/${base_file_name}.all.recal_data.csv
    
    # If we had recalibrated separately by groups, and named each by the
    # group, then we'd had to do the grouping with them instead. For now,
    # since it seems to refuse scattering with group arrays of intervals
    # directly, we are doing all intervals at once and getting thus only
    # a single file that needs not be "gathered".


    # ------------------ ATTENTION -----------------
    # Next we need to apply the BQSR
    BY_SEQ_GRPS='NO'
    if [ "$BY_SEQ_GRPS" == 'YES' ] ; then
	echo "ApplyBQSR by groups NOT YET READY"
        exit
	# DO NOT USE: THE FUNCTION AND LOGIC ARE INCOMPLETE
        cat sequence_grouping_with_unmapped.txt | while read subgroup ; do
            $CALL ApplyBQSR $sorted_all_bam "'$subgroup'"
            recal_all_bam=$align_bwa_dir/$base_file_name.aligned.duplicates_marked.recalibrated.bam
            for sorted_bam in $align_bwa_dir/*.duplicate_marked.sorted.bam ; do
                  $CALL ApplyBQSR $sorted_bam "'$subgroup'"
            done
        done
    else
        echo "Using EXOME.intervals.list + unmapped"
        if [ ! -s sequence_grouping_with_unmapped.intervals.list ] ; then
            cp EXOME.intervals.list sequence_grouping_with_unmapped.intervals.list
	    echo "unmapped" >> sequence_grouping_with_unmapped.intervals.list
        fi
        # Generate the recalibration model by interval
        $CALL ApplyBQSR $sorted_all_bam sequence_grouping_with_unmapped.intervals.list
        for sorted_bam in $align_bwa_dir/*.duplicate_marked.sorted.bam ; do
            $CALL ApplyBQSR $sorted_bam sequence_grouping_with_unmapped.intervals.list
        done
    fi
    
    # Merge the separately corrected BAMS
    # If we work by sequence groups, each SG will generate a separate BAM
    # if we work by individuals with all SG at once, each individual generates
    # a separate BAM
    individual_bams=$( ls $align_bwa_dir/*recalibrated.bam | grep -v $base_file_name )
    
    # Fails with error indexing a read aligned to chr2 (should be 17?)
    ###$CALL GatherBamFiles "${individual_bams[@]}"
    #output_bam=./$align_bwa/$base_file_name.all.align.dup.recal.bam
}


PreProcessingForVariantDiscovery_GATK4
exit

