#!/bin/bash
#
# This script has been developed at Scientific Computing Service of
# the National Biotechnology Center of the Spanish Research Council 
# (CNB-CSIC) and it is therefore the property of CNB-CSIC.
#
# It is based on Cromwell scripts developed by the Broad Institute.
#
# (c) CNB-CSIC. 2019
# Released under a LGPL or EU-LGPL license
#
#---------------------------------------------------------------------
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

set -x
DEBUG='NO'
if [ "$DEBUG" == 'YES' ] ; then
    CALL='echo'
else
    CALL=''
fi

# we should use JAVA-8!!!
export PATH=/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/:$PATH

################### THESE PROVIDE DEFAULT VALUES ###################
################## OVERRIDEN IN gatk4-BP2019.in.sh #################
## SAMPLE NAME AND UNMAPPED BAMS
ref_name=human.genome.19
study_name="patient.data"
flowcell_unmapped_bams_list="./uBam.list"
unmapped_bam_suffix=".bam"
## DATA DIRECTORIES
data_dir="./hg19"
gatk_bundle_dir="$data_dir/gatk-bundle"
## PATHS
bwa_path='/usr/bin/'
gatk_exec="${HOME}/contrib/gatk4/gatk"
gotc_path="${HOME}/contrib/gatk4/"
study_interval_list="$panel_dir/exome.interval_list"
align_dir="./align.gatk"
## "KNOWN SITES RESOURCES"
dbSNP_vcf="$gatk_bundle_dir/dbsnp_138.hg19.vcf"
dbSNP_vcf_index="$gatk_bundle_dir/dbsnp_138.hg19.vcf.idx"
#
known_indels_sites_VCFs=(
    "$gatk_bundle_dir/dbsnp_138.hg19.vcf"
    "$gatk_bundle_dir/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf"
    "$gatk_bundle_dir/1000G_phase1.indels.hg19.sites.vcf"
#	Homo_sapiens_assembly38.known_indels.vcf.gz	N/A for hg19
)

# get actual definitions from a local file
if [ -s gatk4-BP2019.in.sh ] ; then
    . gatk4-BP2019.in.sh
else
    cp `dirname $0`/gatk4-BP2019.in.sh \
       ./gatk4-BP2019.in.sh
    echo "I have created a preferences file in this directory:"
    echo ""
    echo "        gatk4-BP2019.in.sh"
    echo ""
    echo "Please, open it in a text editor, adapt it to your needs and"
    echo "run this command again."
    exit
fi

## REFERENCE FILES
ref_dict="$data_dir/$ref_name.dict"
ref_fasta="$data_dir/$ref_name.fasta"
ref_fasta_index="$data_dir/$ref_name.fasta.fai"
ref_sa="$data_dir/$ref_name.fasta.sa"
ref_amb="$data_dir/$ref_name.fasta.amb"
ref_bwt="$data_dir/$ref_name.fasta.bwt"
ref_ann="$data_dir/$ref_name.fasta.ann"
ref_pac="$data_dir/$ref_name.fasta.pac"

   
## "MISC PROGRAM PARAMETERS", 
compression_level=5
num_cpu="16"


## "JAVA OPTIONS", 
SamToFastqAndBwaMem_java_opt="-Xms3000m"
MergeBamAlignment_java_opt="-Xms3000m"
MarkDuplicates_java_opt="-Xms4000m"
#SortAndFixTags_java_opt_sort="-Xms4000m"
SortAndFixTags_java_opt_sort="-Xms16000m"
#SortAndFixTags_java_opt_fix="-Xms5000m"
SortAndFixTags_java_opt_fix="-Xms16000m"
BQSRecalibrator_java_opt="-Xms4000m"
#GatherBqsrReports_java_opt="-Xms3000m"
GatherBqsrReports_java_opt="-Xms16000m"
ApplyBQSR_java_opt="-Xms3000m"
#GatherBamFiles_java_opt="-Xms2000m"
GatherBamFiles_java_opt="-Xms16000m"
# default
java_opt="-Xms4000m"

## "MEMORY ALLOCATION", 
GetBwaVersion_mem_size="1 GB"
SamToFastqAndBwaMem_mem_size="14 GB"
MergeBamAlignment_mem_size="3500 MB"
MarkDuplicates_mem_size="7 GB"
SortAndFixTags_mem_size="5000 MB"
CreateSequenceGroupingTSV_mem_size="2 GB"
BQSRecalibrator_mem_size="6 GB"
GatherBqsrReports_mem_size="3500 MB"
ApplyBQSR_mem_size="3500 MB"
GatherBamFiles_mem_size="3 GB"
# default
mem_size="7 GB"

## ALL TOGETHER FILE NAME
ensemble_name=$study_name.$ref_name

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


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM 
# for alignment
#
#	we get		original.bam
#	we produce	original.aligned.bam
#
function SamToFastqAndBwaMem () {

    set -o pipefail
    set -e
    
    # input_bam is a full-length absolute path name
    unmapped_bam=$1
    echo ">>> SamToFastqAndBwaMem $unmapped_bam"
    
    # Get the basename, i.e. strip the filepath and the extension
    bam_basename=`basename $unmapped_bam $unmapped_bam_suffix`
    output_unmerged_bam=$align_dir/$bam_basename.aligned.bam
    
    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}

    bwa_commandline="bwa mem -p -v 3 -t ${num_cpu} -Y ${bash_ref_fasta}"

    if [ ! -s ${output_unmerged_bam} ] ; then
	mkdir -p $align_dir
        if [ ! -s $bam_basename.fastq ] ; then
	    echo ">>> Generating $bam_basename.fastq"
            
            # convert to FastQ, pipe to bwa for alignment and to sam for bam comversion
#            java -Dsamjdk.compression_level=${compression_level} \
#    	          ${SamToFastqAndBwaMem_java_opt} \
#                  -jar ${gotc_path}/gatk \
            $gatk_exec \
                  SamToFastq \
	          -INPUT=${unmapped_bam} \
	          -FASTQ=$bam_basename.fastq \
	          -INTERLEAVE=true \
	          -NON_PF=true 
	fi
        echo ">>> Aligning and making ${output_unmerged_bam}"
        ${bwa_path}${bwa_commandline} $bam_basename.fastq   \
            2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) \
            > $bam_basename.sam
        samtools view -1 $bam_basename.sam > ${output_unmerged_bam}
        rm -f $bam_basename.fastq $bam_basename.sam
    fi

}


# Merge original input uBAM file with BWA-aligned BAM file
#
#	we get		original.bam
#	we use		original.aligned.bam
#	we produce	original.aligned.merged.bam
function MergeBamAlignment () {
    unmapped_bam=$1
    
    echo ">>> MergeBamAlignment $unmapped_bam"

    # set the bash variable needed for the command-line
    bash_ref_fasta=${ref_fasta}
    bam_basename=`basename $unmapped_bam $unmapped_bam_suffix`
    unmerged_aligned_bam=$align_dir/$bam_basename.aligned.bam
    output_bam=$align_dir/${bam_basename}.aligned.merged.bam
    
    if [ ! -s $output_bam ] ; then
	echo ">>> Merging BAM alignment with original unmapped bam"
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
# May work with only one or an array of many files.
# If an array is used, then $ensemble_name ($study_name$ref_name) will be 
# used for output combining all of them
#
#	we get		x.aligned.merged.bam
#	we produce	x.aligned.merged.duplicates_marked.bam
function MarkDuplicates() {
    
    input_bams=( $@ )
    
    echo ">>> MarkDuplicates" ${input_bams[@]}
    
    if [ ${#input_bams[@]} -gt 1 ] ; then
        # use a summary name
        #ensemble_name=$study_name.$ref_name	# defined at the top
        output_bam_basename=${ensemble_name}
    else
        # use a file specific name
        input_bam=$input_bams
        output_bam_basename=`basename $input_bam .bam`
    fi
    output_bam=$align_dir/$output_bam_basename.duplicates_marked.bam
    metrics_filename=$align_dir/${output_bam_basename}.duplicate_metrics

    # prepare argument list of input files
    input_list=''
    for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done
    
    
    if [ ! -s ${output_bam} ] ; then
        # This task is assuming query-sorted input so that the Secondary and 
        # Supplementary reads get marked correctly. 
        # This works because the output of BWA is query-grouped and therefore,
        # so is the output of MergeBamAlignment. 
        # While query-grouped isn't actually query-sorted, it's good enough for
        # using MarkDuplicates with ASSUME_SORT_ORDER="queryname"
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${MarkDuplicates_java_opt}" \
            MarkDuplicates \
            ${input_list} \
            --OUTPUT ${output_bam} \
            --METRICS_FILE ${metrics_filename} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true #\
            #--REMOVE_DUPLICATES=true
    fi    
}




# Sort BAM file by coordinate order and fix tag values for NM and UQ
#
#	we get		aligned.merged.duplicates_marked.bam
#	we produce	aligned.merged.duplicate_marked.sorted.bam
#
function SortAndFixTags() {

    input_bam=$1
    
    echo ">>> SortAndFixTags $input_bam"
    
    unsorted_bam_base=`basename $input_bam .bam`
    output_bam=$align_dir/${unsorted_bam_base}.sorted.bam
  
    set -o pipefail

    echo ">>> Sorting $input_bam"
    ### May need to split the command in two, even if it eats up extra space
    if [ ! -s $output_bam ] ; then
    #if [ "NO" == "YES" ] ; then
         ${gatk_exec} \
           --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_sort}" \
           SortSam \
           --TMP_DIR "./scratch" \
           --INPUT ${input_bam} \
           --OUTPUT /dev/stdout \
           --SORT_ORDER "coordinate" \
           --CREATE_INDEX false \
           --CREATE_MD5_FILE false \
         | \
         ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_fix} -Xms24000m" \
           SetNmMdAndUqTags \
           --INPUT /dev/stdin \
           --OUTPUT ${output_bam} \
           --CREATE_INDEX true \
           --CREATE_MD5_FILE true \
           --REFERENCE_SEQUENCE ${ref_fasta}
    fi
    # if the above failed for any reason, try again in two steps
    if [ ! -s $output_bam ] ; then
        if [ ! -s sorted.bam ] ; then
        ${gatk_exec} \
          --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_sort}" \
          SortSam \
          --INPUT ${input_bam} \
          --OUTPUT sorted.bam \
          --SORT_ORDER "coordinate" \
          --CREATE_INDEX false \
          --CREATE_MD5_FILE false 
        fi
        # if that failed, try sorting with SAMTOOLS (which is faster BTW)
        if [ ! -s sorted.bam ] ; then
            echo ">>>" samtools sort -m 8G -o sorted.bam -O bam -T temp.bam -@ 10 ${input_bam}
            samtools sort -m 8G -o sorted.bam -O bam -T temp.bam -@ 10 ${input_bam}
	fi
        if [ ! -s $output_bam ] ; then
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${SortAndFixTags_java_opt_fix} -Xms36000m" \
          SetNmMdAndUqTags \
          --INPUT sorted.bam \
          --OUTPUT ${output_bam} \
          --CREATE_INDEX true \
          --CREATE_MD5_FILE true \
          --REFERENCE_SEQUENCE ${ref_fasta}
        fi
    fi
    echo ">>> Created $output_bam"  ; rm -f sorted.bam
}



# Generate sets of intervals for scatter-gathering over chromosomes
# UNUSED FOR NOW
function CreateSequenceGroupingTSV() {

  if [ ! -s sequence_grouping.txt -o ! -s sequence_grouping_with_unmapped.txt ] 
  then
    # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter. 
    # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
    # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
    ${HOME}/contrib/anaconda2/bin/python <<PYCODE

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

PYCODE
  fi

}

# Generate Base Quality Score Recalibration (BQSR) model
function BQSRecalibrator() {
    input_bam=$1
    shift
    sequence_group_interval=( $@ )

    echo ">>> BQSRecalibrating $input_bam"
    # either we do it for all sequence group intervals at once
    # or we generate one different recalibration_report for each
    # interval, so we can later gather all together in 
    # a call to GatherBsqrReports()
    base_out_name=`basename $input_bam .bam`
    output_recalibration_report=$align_dir/${base_out_name}.recal_data.csv
    
    # known SNPs to consider (using a global variable)
    known_snps_arg="--known-sites ${dbSNP_vcf}"
    
    # Known indels to consider (using a global variable)
    known_indels_arg=''
    for i in ${known_indels_sites_VCFs[@]} ; do
        known_indels_arg="${known_indels_arg} --known-sites $i"
    done
    
    # interval lists from remaining arguments 
    sg_arg=''
    for i in ${sequence_group_interval[@]} ; do
        sg_arg="${sg_arg} -L $i"
    done
    
    if [ ! -s "${output_recalibration_report}" ] ; then
        echo ">>> Creating recalibration reports for $input_bam"
        echo "${gatk_exec} --java-options ${BQSRecalibrator_java_opt} \
          BaseRecalibrator \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${output_recalibration_report} \
          $known_snps_arg \
          $known_indels_arg \
          $sg_arg #\
          #--known-sites ${dbSNP_vcf} \
          #-L sequence_grouping.intervals.list"

        ${gatk_exec} --java-options "${BQSRecalibrator_java_opt}" \
          BaseRecalibrator \
          -R ${ref_fasta} \
          -I ${input_bam} \
          --use-original-qualities \
          -O ${output_recalibration_report} \
          $known_snps_arg \
          $known_indels_arg \
          $sg_arg #\
          #-L sequence_grouping.intervals.list
    fi
}


# Combine multiple recalibration tables from scattered BQSRecalibrator runs
# Note that when run from GATK 3.x the tool is not a walker and is invoked 
# differently.
function GatherBqsrReports {
    input_bqsr_reports=( $@ )
  

    inrep_arg=''
    for i in ${input_bqsr_reports[@]} ; do
        inrep_arg="${inrep_arg} -I $i"
    done
    output_report_filename=$align_dir/${ensemble_name}.all.recal_data.csv

    if [ ! -s ${output_report_filename} ] ; then
        ${gatk_exec} --java-options "${GatherBqsrReports_java_opt}" \
            GatherBQSRReports \
            $inrep_arg \
            -O ${output_report_filename}
    fi
}


# Apply Base Quality Score Recalibration (BQSR) model
# using a previously computed BQSR report (by BQSRecalibrator /
# GatherBqsrReports)

function ApplyBQSR() {
    input_bam=$1
    shift
    sequence_group_interval=( $@ )

    echo ">>> Applying BQSR to $input_bam"
    
    # we use for each file its own specific recalibration report
    base_name=`basename $input_bam .bam`
    recalibration_report=$align_dir/$base_name.recal_data.csv
    output_bam=$align_dir/$base_name.recalibrated.bam

    sg_arg=''
    for i in ${sequence_group_interval[@]} ; do
        sg_arg="${sg_arg} -L $i"
    done
    
    if [ ! -s ${input_bam} ] ; then
    	echo "${input_bam} does not exist!"
        return
    fi
    if [ ! -s ${output_bam} ] ; then
        ${gatk_exec} --java-options "${ApplyBQSR_java_opt}" \
             ApplyBQSR \
             -R ${ref_fasta} \
             -I ${input_bam} \
             -O ${output_bam} \
             $sg_arg \
             -bqsr ${recalibration_report} \
             --static-quantized-quals 10 \
	     --static-quantized-quals 20 \
	     --static-quantized-quals 30 \
             --add-output-sam-program-record \
             --create-output-bam-md5 \
             --use-original-qualities	# what is the sense of this last arg?
    fi
}


# Combine multiple recalibrated BAM files from scattered ApplyRecalibration 
# runs
# This should also allow to gather all duplicates_marked.bam files into a
# single one instead of doing MarkDuplicates for all
function GatherBamFiles() {
    input_bams=( $@ )
  
    in_bams_arg=''
    for i in ${input_bams[@]} ; do  
        in_bams_arg="${in_bams_arg} --INPUT $i"
    done
    output_bam=./$align_dir/$ensemble_name.all.recal.bam

    if [ ! -s $output_bam ] ; then
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${GatherBamFiles_java_opt}" \
             GatherBamFiles \
             ${in_bams_arg} \
             --OUTPUT ${output_bam} \
             --CREATE_INDEX true \
             --CREATE_MD5_FILE true
    fi
}


############################################################################
#+------------------------------------------------------------------------+#
#|                                                                        |#
#|                        WORKFLOW DEFINITION                             |#
#|                                                                        |#
#+------------------------------------------------------------------------+#
############################################################################


function PreProcessingForVariantDiscovery_GATK4() {
    # Get the version of BWA to include in the PG record in the header of the BAM produced 
    # by MergeBamAlignment. 
    bwa_version=`GetBwaVersion`

    # Align flowcell-level unmapped input bams
    #for unmapped_bam in "${flowcell_unmapped_bams[@]}" ; do
    for unmapped_bam in $(cat $flowcell_unmapped_bams_list) ; do

	echo ">>> Processing $unmapped_bam"

        # Map reads to reference
        $CALL SamToFastqAndBwaMem $unmapped_bam
        # generates $unmapped.aligned.bam

        # Merge original uBAM and BWA-aligned BAM 
        $CALL MergeBamAlignment $unmapped_bam
        # generates $unmapped.aligned.merged.bam
        
    done
    
    # in former BP protocols, after BWA, the sam file was sorted first
    # and the deduplicated. Here we remove duplicates first and sort later.
    
    # get list of bam files generated into an array
    # we have generated a file "X.aligned.bam" and an "X.aligned.mergedb.am"
    # in $align_dir
    merged_bams=$( ls $align_dir/*.merged.bam ) 
    
    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    #MergeBamALignment $merged_bams # not needed we'll do it with MarkDuplicates

    # mark duplicates on a per-file basis
    for i in ${merged_bams[@]} ; do
        out=$align_dir/`basename $i .bam`.aligned.unsorted.duplicates_marked.bam
        echo ">>> Marking Duplicates in $out"
        if [ ! -s $out ] ; then
            $CALL MarkDuplicates $i
	fi
    done

    # We take advantage of the tool's ability to take multiple BAM inputs and 
    # write out a single output to avoid having to spend time just 
    # merging BAM files.
    ### JR ### This call should do it for all at once and generate
    # a single file, $ensemble_name.duplicates_marked.bam
    ### note that it may fail for many very big files, in that case, we
    ### probably should consider using MergeBamAlignment
    $CALL MarkDuplicates "${merged_bams[@]}"


    # Sort deduplicated BAM files and fix tags
    # This will also consider the aggregate file if it exists
    for unsorted in $align_dir/*.duplicates_marked.bam ; do
        echo "I'm going to sort $unsorted"
        sorted_bam=$align_dir/`basename $unsorted .bam`.sorted.bam
        if [ ! -e $sorted_bam ] ; then
            $CALL SortAndFixTags $unsorted
        fi
        # this generates "$unsorted.sorted.bam"
        # it will also sort the ensemble file if it exists
    done


    #if [ "YES" == "WGS" ] ; then
    #     echo -n ""
    #    ### This was used in the original workflow to split and reunite in parts a
    #    ### large BAM file for parallel processing.
    #    ### Assuming that all samples were in a single, huge file.
    #
    #    # Create list of sequences for scatter-gather parallelization 
    #    CreateSequenceGroupingTSV 
    #    # read groupings into arrays
    #    sequence_grouping=$(< "sequence_grouping.txt")
    #    sequence_grouping_with_unmapped=$(< "sequence_grouping_with_unmapped.txt")
    #
    #    # Here we may want to use our own groupings only for the genes we are
    #    # interested in, saving the global groupings and generating new ones:
    #    # We will group by chromosome
    #    if [ ! -e EXOME.intervals.list ] ; then
    #        echo "Please, provide an EXOME.intervals.list file"
    #        echo "e.g. ln -s ./$panel_dir/INGEM_panel5_120_mod_interval_format.intervals EXOME.intervals.list"
    #        exit
    #    fi
    #    USE_EXOME_COORDINATES='yes'
    #    if [ "$USE_EXOME_COORDINATES" == "yes" ] ; then
    #        if [ ! -e sequence_grouping.txt.sav ] ; then
    #            mv sequence_grouping.txt sequence_grouping.txt.sav
    #            mv sequence_grouping_with_unmapped.txt sequence_grouping_with_unmapped.txt.sav
    #            echo -n "" > sequence_grouping.txt
    #            cut -d: -f1 EXOME.intervals.list | sort | uniq |\
    #                while read chr ; do
    #                    echo "grouping $chr"
    #                    grep $chr EXOME.intervals.list | \
    #                    while IFS=$'\n' read interval ; do
    #                        echo -n "$interval	" >> sequence_grouping.txt
    #                    done
    #                    echo "">> sequence_grouping.txt
    #                done
    #            cp sequence_grouping.txt sequence_grouping_with_unmapped.txt
    #            echo "unmapped" >> sequence_grouping_with_unmapped.txt
    #        fi
    #    fi
    #fi

    # ------------------ ATTENTION -----------------    
    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM
    
    echo ">>> BQSR Using $study_interval_list"
    # Generate the recalibration model by interval
    # this will also include the ensemble
    for sorted_bam in $align_dir/*.sorted.bam ; do
        $CALL BQSRecalibrator $sorted_bam $study_interval_list
    done
    # get report name for full report (all samples together)
    ### JR ### Disable to disable all together
#    base_report_name=`basename $sorted_all_bam .aligned.duplicate_marked.sorted.bam`
#    recalibration_report=$align_dir/${base_report_name}.recal_data.csv

    # Merge the individual recalibration reports resulting from 
    # recalibration.
    # If we recalibrated by sequence groups, each one will produce a
    # recalibration report. If we did by individuals, recalibrating
    # each with all SGs at once, then each individual will produce a
    # recalibration report. THIS IS OUR CASE
    # We have
    #	- a recalibration report for all individuals merged with all SG
    #   - a report for each individual with all SG at once
    #	We'll gather the later individual reports into one report
    individual_reports=$( ls $align_dir/*.recal_data.csv | grep -v $ensemble_name )
    
    ### JR ### Disable to disable all together
    $CALL GatherBqsrReports "${individual_reports[@]}"
    ensemble_bqsr_report=$align_dir/${ensemble_name}.all.recal_data.csv
    

    # ------------------ ATTENTION -----------------
    # Next we need to apply the BQSR
    echo ">>> Using $study_intervals + unmapped"
    # the next lines are mostly mnemotechnic sugar
    study_intervals_name=${study_intervals%.*}	# delete shortest match from back
    study_intervals_type=${study_intervals##*.}	# delete longest match from front
    study_useless_edit=$study_intervals/old/new} 	# substitute first old by new 
    study_useless_edit=$study_intervals//old/new} 	# substitute all old by new 
    sequence_groups_plus_unmapped_intervals=${study_intervals_name}+unmapped.${study_intervals_type}

    if [ ! -s $sequence_groups_plus_unmapped_intervals ] ; then
        cp $study_intervals $sequence_groups_plus_unmapped_intervals
	echo "unmapped" >> $sequence_groups_plus_unmapped_intervals
    fi
    # Generate the recalibration model by interval for all sorted files
    # this includes the ensemble file
    for sorted_bam in $align_dir/*.sorted.bam ; do
        $CALL ApplyBQSR $sorted_bam $sequence_groups_plus_unmapped_intervals
    done

    
    # Merge the separately corrected BAMS
    # If we work by sequence groups, each SG will generate a separate BAM
    # if we work by individuals with all SG at once, each individual 
    # corresponds to a separate BAM
    #individual_bams=$( ls $align_dir/*.recalibrated.bam | grep -v $ensemble_name )
    
    # This should generate a new BAM equivalent of the ensemble BAM.
    # Fails with error indexing a read aligned to chr2 (should be 17?)
    ###$CALL GatherBamFiles "${individual_bams[@]}"
    #output_bam=./$align_dir/$ensemble_name.all.align.dup.recal.bam
}


PreProcessingForVariantDiscovery_GATK4
exit

