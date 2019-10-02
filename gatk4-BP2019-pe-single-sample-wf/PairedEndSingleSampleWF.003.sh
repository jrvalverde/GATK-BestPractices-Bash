ref=ucsc.hg19
ref_fasta=$ref.fasta
ref_fasta_index=$ref_fasta.fai
ref_dict=$ref.dict
# This is the .alt file from bwa-kit (https://github.com/lh3/bwa/tree/master/bwakit),
# listing the reference contigs that are "alternative".
ref_alt=""####<<<<####


tools=/usr/bin
gatk3=$HOME/contrib/gatk3
gatk4=$HOME/contrib/gatk4
gatk_exec=$gatk4/gatk
#picard_exec="java -jar $gatk3/picard.jar"
picard_exec=${gatk_exec}
VerifyBamID_exec=$gatk4/contrib/verifyBamID

  # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
  # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
  max_duplication_in_reasonable_sample=0.30
  max_chimerism_in_reasonable_sample=0.15
  bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
  compression_level=2	# for SAM jdk



function PairedEndSingleSampleWorkflow() {

  #global File contamination_sites_ud
  #global File contamination_sites_bed
  #global File contamination_sites_mu
  #global File? fingerprint_genotypes_file
  #global File? fingerprint_genotypes_index
  #global File? haplotype_database_file
  #global File wgs_evaluation_interval_list
  #global File wgs_coverage_interval_list

  #String sample_name
  #String base_file_name
  #String final_gvcf_base_name
  #File flowcell_unmapped_bams_list
  #Array[File] flowcell_unmapped_bams = read_lines(flowcell_unmapped_bams_list)
  #String unmapped_bam_suffix

  #File wgs_calling_interval_list
  #Int haplotype_scatter_count
  #Int break_bands_at_multiples_of

  #Int? read_length

  #global File ref_fasta
  #global File ref_fasta_index
  #global File ref_dict
  #global File ref_alt
  #global File ref_bwt
  #global File ref_sa
  #global File ref_amb
  #global File ref_ann
  #global File ref_pac

  #global File dbSNP_vcf
  #global File dbSNP_vcf_index
  #global Array[File] known_indels_sites_VCFs
  #global Array[File] known_indels_sites_indices

    ### Provide some defaults in case the variables are undefined for
    ### whatever reason.

    # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
    # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
    #Float max_duplication_in_reasonable_sample = 0.30
    #Float max_chimerism_in_reasonable_sample = 0.15
    max_duplication_in_reasonable_sample=${max_duplication_in_reasonable_sample:-"0.30"}
    max_chimerism_in_reasonable_sample=${max_chimerism_in_reasonable_sample:-"0.15"}

    #String bwa_commandline="bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"
    bwa_commandline=${bwa_commandline:-"bwa mem -K 100000000 -p -v 3 -t 16 -Y $bash_ref_fasta"}

    #global Int compression_level = 2
    compression_level=${compression_level:-2}


    ### Start processing

    # Get the version of BWA to include in the PG record in the header of the BAM produced
    # by MergeBamAlignment.
    bwa_version=`GetBwaVersion`

    # Align flowcell-level unmapped input bams 

    # QC the unmapped BAM
    CollectQualityYieldMetrics $unmapped_bam
    #  input:
    #    input_bam = unmapped_bam,
    #    metrics_filename = unmapped + ".unmapped.quality_yield_metrics"
    yield_metrics=${unmapped_bam/%.bam/.unmapped.quality_yield_metrics}


    # Map reads to reference
    SamToFastqAndBwaMemAndMba $unmapped_bam
    #  input:
    #    input_bam = unmapped_bam,
    #    output_bam_basename = unmapped + ".aligned",
    #    global bwa_commandline = bwa_commandline,
    #    global ref_fasta = ref_fasta,
    #    global ref_fasta_index = ref_fasta_index,
    #    global ref_dict = ref_dict,
    #    global ref_alt = ref_alt,
    #    global ref_bwt = ref_bwt,
    #    global ref_amb = ref_amb,
    #    global ref_ann = ref_ann,
    #    global ref_pac = ref_pac,
    #    global ref_sa = ref_sa,
    #    bwa_version = GetBwaVersion.version,
    #    global compression_level = compression_level,
    aligned_bam=${unmapped_bam/%.bam/.aligned.bam}


    # QC the aligned but unsorted readgroup BAM
    # we use no reference as the input here is unsorted, providing a reference 
    # would cause an error
    CollectUnsortedReadgroupBamQualityMetrics $aligned_bam 
    #  input:
    #    input_bam = SamToFastqAndBwaMemAndMba.output_bam,
    #    output_bam_prefix = aligned + ".readgroup",
    dist_metrics=${aligned_bam/%.bam/.readgroup}.quality_distribution_metrics


    # Aggregate aligned+merged flowcell BAM files and mark duplicates
    # We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
    # to avoid having to spend time just merging BAM files.
    MarkDuplicates $aligned_bam
    #  input:
    #    input_bams = SamToFastqAndBwaMemAndMba.output_bam,
    #    output_bam_basename = base_file_name + ".aligned.unsorted.duplicates_marked",
    #    metrics_filename = base_file_name + ".duplicate_metrics",
    #    global compression_level = compression_level,
    marked_bam=${aligned_bam%.bam}.duplicates_marked.bam
    duplicate_metrics=${aligned_bam%.bam}.duplicate_metrics


    # Sort aggregated+deduped BAM file and fix tags
    SortSam $marked_bam
    #input:
    #  input_bam = MarkDuplicates.output_bam,
    #  output_bam_basename = base_file_name + ".aligned.duplicate_marked.sorted",
    #  global compression_level = compression_level,
    sorted_bam=${marked_bam/%.bam/.sorted.bam}

    if [ "$haplotype_database_file" != "" ] then
        # Check identity of fingerprints across readgroups
        CrossCheckFingerprints $sorted_bam
        #  input:
        #    input_bams = SortSampleBam.output_bam,
        #    input_bam_indexes = SortSampleBam.output_bam_index,
        #    global haplotype_database_file = haplotype_database_file,
        #    metrics_filename = base_file_name + ".crosscheck"
        ccfp_metrics=${sorted_bam/%.bam/.crosscheck}
    fi

    # Create list of sequences for scatter-gather parallelization
    CreateSequenceGroupingTSV
    #  input:
    #    ref_dict = ref_dict
    sequence_grouping="sequence_grouping.txt"
    sequence_grouping_with_unmapped="sequence_grouping_with_unmapped.txt"

    # Estimate level of cross-sample contamination
    contamination=`CheckContamination $sorted_bam 0.75`
    #input:
    #  input_bam = SortSampleBam.output_bam,
    #  input_bam_index = SortSampleBam.output_bam_index,
    #  global contamination_sites_ud = contamination_sites_ud,
    #  global contamination_sites_bed = contamination_sites_bed,
    #  global contamination_sites_mu = contamination_sites_mu,
    #  global ref_fasta = ref_fasta,
    #  global ref_fasta_index = ref_fasta_index,
    #  output_prefix = base_file_name + ".preBqsr",
    #  contamination_underestimation_factor = 0.75
    # NOTE: it will also produce some output on stdout
    selfSM=${sorted_bam%.bam}.selfSM

    # Perform Base Quality Score Recalibration (BQSR) on the sorted BAM
#  NOT SCATTERED, we'll do it all at once 
#  scatter (subgroup in CreateSequenceGroupingTSV.sequence_grouping) {
    #
    # Generate the recalibration model
    intervals=$study_interval_list
    # XXX JR XXX NOTE: may require reading the file into an array as
    #    intervals=($(<$study_interval_list))
    BaseRecalibrator $sorted_bam $intervals
    #  input:
    #    input_bam = SortSampleBam.output_bam,
    #    recalibration_report_filename = base_file_name + ".recal_data.csv",
    #    ###sequence_group_interval = subgroup,
    #    global dbSNP_vcf = dbSNP_vcf,
    #    global dbSNP_vcf_index = dbSNP_vcf_index,
    #    global known_indels_sites_VCFs = known_indels_sites_VCFs,
    #    global known_indels_sites_indices = known_indels_sites_indices,
    #    global ref_dict = ref_dict,
    #    global ref_fasta = ref_fasta,
    #    global ref_fasta_index = ref_fasta_index,
    recalibration_report=${sorted_bam%.bam}.recal_data.csv

#  # NOT DONE AS WE ARE GENERATING THE REPORT ALL AT ONCE
#  # Merge the recalibration reports resulting from by-interval recalibration
#  # The reports are always the same size
#  call GatherBqsrReports {
#    input:
#      input_bqsr_reports = BaseRecalibrator.recalibration_report,
#      output_report_filename = base_file_name + ".recal_data.csv",

    # Apply the recalibration model
    ApplyBQSR $sorted_bam $recalibration_report $intervals 
    #  input:
    #    input_bam = SortSampleBam.output_bam,
    #    output_bam_basename = recalibrated_bam_basename,
    #    recalibration_report = GatherBqsrReports.output_bqsr_report,
    #    sequence_group_interval = intervals,
    #    global ref_dict = ref_dict,
    #    global ref_fasta = ref_fasta,
    #    global ref_fasta_index = ref_fasta_index,
    #    global compression_level = compression_level,
    recal_bam=${input_bam%.bam}.recalibrated.bam
    recal_report= ${input_bam%.bam}
  
  
#  # UNNEEDED, we did it all at once
#  # Merge the recalibrated BAM files resulting from by-interval recalibration
#  call GatherBamFiles {
#    input:
#      input_bams = ApplyBQSR.recalibrated_bam,
#      output_bam_basename = base_file_name,
#      # Multiply the input bam size by two to account for the input and output
#      disk_size = (2 * agg_bam_size) + additional_disk,
#      compression_level = compression_level,
#      preemptible_tries = agg_preemptible_tries,
#      docker = gitc_docker,
#      gitc_path = gitc_path
#  }
#


    # QC the final BAM (consolidated after scattered BQSR)
    CollectReadgroupBamQualityMetrics $recal_bam
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  input_bam_index = GatherBamFiles.output_bam_index,
    #  output_bam_prefix = base_file_name + ".readgroup",
    #  global ref_dict = ref_dict,
    #  global ref_fasta = ref_fasta,
    #  global ref_fasta_index = ref_fasta_index,


    # QC the final BAM some more (no such thing as too much QC)
    CollectAggregationMetrics $recal_bam
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  input_bam_index = GatherBamFiles.output_bam_index,
    #  output_bam_prefix = base_file_name,
    #  global ref_dict = ref_dict,
    #  global ref_fasta = ref_fasta,
    #  global ref_fasta_index = ref_fasta_index,
    alignment_summary=${recal_bam%.bam}.alignment_summary_metrics

    if [ "$haplotype_database_file" != "" -a "$fingerprint_genotypes_file" != "" ]
    then
        # Check the sample BAM fingerprint against the sample array
        sample=$sample_name
        CheckFingerprint $recal_bam $sample
        #  input:
        #    input_bam = GatherBamFiles.output_bam,
        #    input_bam_index = GatherBamFiles.output_bam_index,
        #    global haplotype_database_file = haplotype_database_file,
        #    global genotypes = fingerprint_genotypes_file,
        #    global genotypes_index = fingerprint_genotypes_index,
        #    output_basename = base_file_name,
        #    sample = sample_name,
        fingerprint_summary=${recal_bam%.bam}.fingerprinting_summary_metrics
    fi

    # QC the sample WGS metrics (stringent thresholds)
    # NOTE: we'll do this check even if when are working with exome data
    read_size=${read_length:-250}
    interval_list=$study_interval_list
    CollectWgsMetrics $recal.bam $interval_list $read_size
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  input_bam_index = GatherBamFiles.output_bam_index,
    #  metrics_filename = base_file_name + ".wgs_metrics",
    #  ref_fasta = ref_fasta,
    #  ref_fasta_index = ref_fasta_index,
    #  wgs_coverage_interval_list = wgs_coverage_interval_list,
    #  read_length = read_length,
    wgs_metrics=${recal_bam%.bam}.wgs_metrics

    # QC the sample raw WGS metrics (common thresholds)
    CollectRawWgsMetrics $recal_bam $interval_list $read_size
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  input_bam_index = GatherBamFiles.output_bam_index,
    #  metrics_filename = base_file_name + ".raw_wgs_metrics",
    #  ref_fasta = ref_fasta,
    #  ref_fasta_index = ref_fasta_index,
    #  wgs_coverage_interval_list = wgs_coverage_interval_list,
    #  read_length = read_length,
    raw_wgs_metrics=${recal_bam%.bam}.raw_wgs_metrics
    
    
    # Generate a checksum per readgroup in the final BAM
    CalculateReadGroupChecksum $recal_bam
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  input_bam_index = GatherBamFiles.output_bam_index,
    #  read_group_md5_filename = recalibrated_bam_basename + ".bam.read_group_md5",
    recal_md5=${recal_bam}.read_group_md5

    # Convert the final merged recalibrated BAM file to CRAM format
    ConvertToCram $recal_bam
    #input:
    #  input_bam = GatherBamFiles.output_bam,
    #  global ref_fasta = ref_fasta,
    #  global ref_fasta_index = ref_fasta_index,
    #  output_basename = base_file_name,
    recal_cram=${recal_bam%.bam}.cram

    # ValidateSamFile runs out of memory in mate validation on crazy edge case data, so we want to skip the mate validation
    # in those cases.  These values set the thresholds for what is considered outside the normal realm of "reasonable" data.
    #Float max_duplication_in_reasonable_sample = 0.30
    #Float max_chimerism_in_reasonable_sample = 0.15
    # These values have already been provided as default above
    max_duplication_in_reasonable_sample=0.30
    max_chimerism_in_reasonable_sample=0.15

    # Check whether the data has massively high duplication or chimerism rates
    CheckPreValidation $duplicate_metrics \
        $alignment_summary \
  	$max_duplication_in_reasonable_sample \
        $max_chimerism_in_reasonable_sample
    #input:
    #  duplication_metrics = MarkDuplicates.duplicate_metrics,
    #  chimerism_metrics = CollectAggregationMetrics.alignment_summary_metrics,
    #  max_duplication_in_reasonable_sample = max_duplication_in_reasonable_sample,
    #  max_chimerism_in_reasonable_sample = max_chimerism_in_reasonable_sample,
    duplication_rate=`cat duplication_value.txt`
    chimerism_rate=`cat chimerism_value.txt`
    if (( $(echo "$duplication_rate > $max_duplication_in_reasonable_sample || $chimerism_rate > $max_chimerism_in_reasonable_sample" |bc -l) )); then
        is_outlier_data="true"
    else
        is_outlier_data="false"
    fi

    # Validate the CRAM file
    ignore="MISSING_TAG_NM"
    max_output=1000000000
    ValidateSamFile $recal_cram $is_outlier_data $max_output $ignore 
    #   input:
    #    input_bam = ConvertToCram.output_cram,
    #    input_bam_index = ConvertToCram.output_cram_index,
    #    report_filename = base_file_name + ".cram.validation_report",
    #    global ref_dict = ref_dict,
    #    global ref_fasta = ref_fasta,
    #    global ref_fasta_index = ref_fasta_index,
    #    ignore = ["MISSING_TAG_NM"],
    #    max_output = 1000000000,
    #    is_outlier_data = CheckPreValidation.is_outlier_data,
    validation_report=${recal_cram}.validation_report
    

    # Break the calling interval_list into sub-intervals
    # Perform variant calling on the sub-intervals, and then gather the results
    ScatterIntervalList $wgs_calling_interval_list \
    			$haplotype_scatter_count \
                        $break_bands_at_multiples_of
    #  input:
    #    interval_list = wgs_calling_interval_list,
    #    scatter_count = haplotype_scatter_count,
    #    break_bands_at_multiples_of = break_bands_at_multiples_of,


  # We need disk to localize the sharded input and output due to the scatter for HaplotypeCaller.
  # If we take the number we are scattering by and reduce by 20 we will have enough disk space
  # to account for the fact that the data is quite uneven across the shards.
  #Int potential_hc_divisor = ScatterIntervalList.interval_count - 20
  #Int hc_divisor = if potential_hc_divisor > 1 then potential_hc_divisor else 1

  # Call variants in parallel over WGS calling intervals
#  scatter (index in range(ScatterIntervalList.interval_count)) {
    # Generate GVCF by interval
    # NOTE: what we do here
    #	On an exome dataset, we use an interval list that contains
    #	all the sequenced intervals
    #	If this were WGS, we could generate an interval list above and
    #   then read it on an array, and use a for loop to run HC on each
    #   of the generated intervals, or just use an IL with all the chromosomes
    #	and use it to do all at once.
    interval_list=$study_interval_list    

    HaplotypeCaller $recal_bam $interval_list $contamination
    #  input:
    #    contamination = CheckContamination.contamination,
    #    input_bam = GatherBamFiles.output_bam,
    #    interval_list = ScatterIntervalList.out[index],
    #    local gvcf_basename = base_file_name,
    #    global ref_dict = ref_dict,
    #    global ref_fasta = ref_fasta,
    #    global ref_fasta_index = ref_fasta_index,
    gvcf=${recal_bam%.bam}.vcf.gz

#  #UNNEEDED we did all at once
#  # Combine by-interval GVCFs into a single sample GVCF file
#  call MergeVCFs {
#    input:
#      input_vcfs = HaplotypeCaller.output_gvcf,
#      input_vcfs_indexes = HaplotypeCaller.output_gvcf_index,
#      output_vcf_name = final_gvcf_base_name + ".g.vcf.gz",


    # Validate the GVCF output of HaplotypeCaller
    ValidateGVCF $gvcf $study_interval_list
    #input:
    #  input_vcf = MergeVCFs.output_vcf,
    #  input_vcf_index = MergeVCFs.output_vcf_index,
    #  global dbSNP_vcf = dbSNP_vcf,
    #  global dbSNP_vcf_index = dbSNP_vcf_index,
    #  global ref_fasta = ref_fasta,
    #  global ref_fasta_index = ref_fasta_index,
    #  global ref_dict = ref_dict,
    #  global wgs_calling_interval_list = wgs_calling_interval_list,


    # QC the GVCF
    CollectGvcfCallingMetrics $gvcf $study_interval_list
    #  input:
    #    input_vcf = MergeVCFs.output_vcf,
    #    input_vcf_index = MergeVCFs.output_vcf_index,
    #    metrics_basename = final_gvcf_base_name,
    #    global dbSNP_vcf = dbSNP_vcf,
    #    global dbSNP_vcf_index = dbSNP_vcf_index,
    #    global ref_dict = ref_dict,
    #    wgs_evaluation_interval_list = wgs_evaluation_interval_list,

  ### LIST OF OUTPUTS
  # Outputs that will be retained when execution is complete
  #  output {
  #    Array[File] quality_yield_metrics = CollectQualityYieldMetrics.metrics
  #
  #    Array[File] unsorted_read_group_base_distribution_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_pdf
  #    Array[File] unsorted_read_group_base_distribution_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.base_distribution_by_cycle_metrics
  #    Array[File] unsorted_read_group_insert_size_histogram_pdf = CollectUnsortedReadgroupBamQualityMetrics.insert_size_histogram_pdf
  #    Array[File] unsorted_read_group_insert_size_metrics = CollectUnsortedReadgroupBamQualityMetrics.insert_size_metrics
  #    Array[File] unsorted_read_group_quality_by_cycle_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_pdf
  #    Array[File] unsorted_read_group_quality_by_cycle_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_by_cycle_metrics
  #    Array[File] unsorted_read_group_quality_distribution_pdf = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_pdf
  #    Array[File] unsorted_read_group_quality_distribution_metrics = CollectUnsortedReadgroupBamQualityMetrics.quality_distribution_metrics
  #
  #    File read_group_alignment_summary_metrics = CollectReadgroupBamQualityMetrics.alignment_summary_metrics
  #    File read_group_gc_bias_detail_metrics = CollectReadgroupBamQualityMetrics.gc_bias_detail_metrics
  #    File read_group_gc_bias_pdf = CollectReadgroupBamQualityMetrics.gc_bias_pdf
  #    File read_group_gc_bias_summary_metrics = CollectReadgroupBamQualityMetrics.gc_bias_summary_metrics
  #
  #    File? cross_check_fingerprints_metrics = CrossCheckFingerprints.metrics
  #
  #    File selfSM = CheckContamination.selfSM
  #    Float contamination = CheckContamination.contamination
  #
  #    File calculate_read_group_checksum_md5 = CalculateReadGroupChecksum.md5_file
  #
  #    File agg_alignment_summary_metrics = CollectAggregationMetrics.alignment_summary_metrics
  #    File agg_bait_bias_detail_metrics = CollectAggregationMetrics.bait_bias_detail_metrics
  #    File agg_bait_bias_summary_metrics = CollectAggregationMetrics.bait_bias_summary_metrics
  #    File agg_gc_bias_detail_metrics = CollectAggregationMetrics.gc_bias_detail_metrics
  #    File agg_gc_bias_pdf = CollectAggregationMetrics.gc_bias_pdf
  #    File agg_gc_bias_summary_metrics = CollectAggregationMetrics.gc_bias_summary_metrics
  #    File agg_insert_size_histogram_pdf = CollectAggregationMetrics.insert_size_histogram_pdf
  #    File agg_insert_size_metrics = CollectAggregationMetrics.insert_size_metrics
  #    File agg_pre_adapter_detail_metrics = CollectAggregationMetrics.pre_adapter_detail_metrics
  #    File agg_pre_adapter_summary_metrics = CollectAggregationMetrics.pre_adapter_summary_metrics
  #    File agg_quality_distribution_pdf = CollectAggregationMetrics.quality_distribution_pdf
  #    File agg_quality_distribution_metrics = CollectAggregationMetrics.quality_distribution_metrics
  #
  #    File? fingerprint_summary_metrics = CheckFingerprint.summary_metrics
  #    File? fingerprint_detail_metrics = CheckFingerprint.detail_metrics
  #
  #    File wgs_metrics = CollectWgsMetrics.metrics
  #    File raw_wgs_metrics = CollectRawWgsMetrics.metrics
  #
  #    File gvcf_summary_metrics = CollectGvcfCallingMetrics.summary_metrics
  #    File gvcf_detail_metrics = CollectGvcfCallingMetrics.detail_metrics
  #
  #    File duplicate_metrics = MarkDuplicates.duplicate_metrics
  #    File output_bqsr_reports = GatherBqsrReports.output_bqsr_report
  #
  #    File output_cram = ConvertToCram.output_cram
  #    File output_cram_index = ConvertToCram.output_cram_index
  #    File output_cram_md5 = ConvertToCram.output_cram_md5
  #
  #    File validate_cram_file_report = ValidateCram.report
  #
  #    File output_vcf = MergeVCFs.output_vcf
  #    File output_vcf_index = MergeVCFs.output_vcf_index
  #  }
}









# Collect sequencing yield quality metrics
function CollectQualityYieldMetrics() {
  if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  
  unmapped_bam=$1		# this should be the unmapped BAM file

  # substitute trailing .bam by .metrics
  metrics_filename=${unmapped_bam/%.bam/.unmapped.quality_yield_metrics}
  
  if [ ! -s $metrics_filename ] ; then {
    #java -Xms2000m -jar $gatk3/picard.jar \
    ${picard_exec} -Xms2000m \
      CollectQualityYieldMetrics \
      INPUT=${unmapped_bam} \
      OQ=true \
      OUTPUT=${metrics_filename}
  }
  metrics="${metrics_filename}"
}


# Get version of BWA
function GetBwaVersion() {
    $tools/bwa 2>&1 | \
    grep -e '^Version' | \
    sed 's/Version: //'
    # Output will be printed to (and collected from) stdout!!!
    # That is why here we do not produce any other output.
}


# Read unmapped BAM, convert on-the-fly to FASTQ and stream to BWA MEM for alignment, then stream to MergeBamAlignment
function SamToFastqAndBwaMemAndMba() {
    if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
    else echo "${FUNCNAME[0]}" $* ; fi

    unmapped_bam=$1	# this should be the unmapped BAM file

    bwa_version=`GetBwaVersion`
    output_bam_basename=${unmapped_bam/%.bam/.aligned}

    # if ref_alt has data in it,
    if [ ! -s ${ref_alt} ]; then exit 1 ; fi

    (
      set -o pipefail
      set -e

      # set the bash variable needed for the command-line
      bash_ref_fasta=${ref_fasta}
      #java -Xms5000m -jar $gatk3/picard.jar \
      ${picard_exec} \
        --java-options "-Xms5000m" \
        SamToFastq \
        INPUT=${unmapped_bam} \
        FASTQ=/dev/stdout \
        INTERLEAVE=true \
        NON_PF=true | \
      ${tools}/${bwa_commandline} /dev/stdin - 2> >(tee ${output_bam_basename}.bwa.stderr.log >&2) | \
      ${picard_exec} \
        --java-options "-Dsamjdk.compression_level=${compression_level} \
          -Xms3000m" \
        MergeBamAlignment \
        VALIDATION_STRINGENCY=SILENT \
        EXPECTED_ORIENTATIONS=FR \
        ATTRIBUTES_TO_RETAIN=X0 \
        ATTRIBUTES_TO_REMOVE=NM \
        ATTRIBUTES_TO_REMOVE=MD \
        ALIGNED_BAM=/dev/stdin \
        UNMAPPED_BAM=${unmapped_bam} \
        OUTPUT=${output_bam_basename}.bam \
        REFERENCE_SEQUENCE=${ref_fasta} \
        PAIRED_RUN=true \
        SORT_ORDER="unsorted" \
        IS_BISULFITE_SEQUENCE=false \
        ALIGNED_READS_ONLY=false \
        CLIP_ADAPTERS=false \
        MAX_RECORDS_IN_RAM=2000000 \
        ADD_MATE_CIGAR=true \
        MAX_INSERTIONS_OR_DELETIONS=-1 \
        PRIMARY_ALIGNMENT_STRATEGY=MostDistant \
        PROGRAM_RECORD_ID="bwamem" \
        PROGRAM_GROUP_VERSION="${bwa_version}" \
        PROGRAM_GROUP_COMMAND_LINE="${bwa_commandline}" \
        PROGRAM_GROUP_NAME="bwamem" \
        UNMAPPED_READ_STRATEGY=COPY_TO_TAG \
        ALIGNER_PROPER_PAIR_FLAGS=true \
        UNMAP_CONTAMINANT_READS=true \
        ADD_PG_TAG_TO_READS=false

      grep -m1 "read .* ALT contigs" ${output_bam_basename}.bwa.stderr.log | \
      grep -v "read 0 ALT contigs"

    )
    # Output:
    # output_bam_basename=unmapped.aligned
    # File output_bam = "${output_bam_basename}.bam"
    # File bwa_stderr_log = "${output_bam_basename}.bwa.stderr.log"
    # the final grep output will go to stdout!!!
  }
}


# Collect base quality and insert size metrics
# QC the aligned but unsorted readgroup BAM
# no reference as the input here is unsorted, providing a reference would cause an error
function CollectUnsortedReadgroupBamQualityMetrics() {
  if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  
  aligned_bam=$1		# this should be the aligned.bam file

  # substitute trailing .bam by .readgroup
  output_bam_prefix=${aligned_bam/%.bam/.readgroup}

  if [ ! -s ${output_bam_prefix}.quality_distribution_metrics ] then
    ${picard_exec} \
      --java-options "-Xms5000m" \
      CollectMultipleMetrics \
      INPUT=${aligned_bam} \
      OUTPUT=${output_bam_prefix} \
      ASSUME_SORTED=true \
      PROGRAM="null" \
      PROGRAM="CollectBaseDistributionByCycle" \
      PROGRAM="CollectInsertSizeMetrics" \
      PROGRAM="MeanQualityByCycle" \
      PROGRAM="QualityScoreDistribution" \
      METRIC_ACCUMULATION_LEVEL="null" \
      METRIC_ACCUMULATION_LEVEL="ALL_READS"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf
  fi

#  output {
#    File base_distribution_by_cycle_pdf = "${output_bam_prefix}.base_distribution_by_cycle.pdf"
#    File base_distribution_by_cycle_metrics = "${output_bam_prefix}.base_distribution_by_cycle_metrics"
#    File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
#    File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
#    File quality_by_cycle_pdf = "${output_bam_prefix}.quality_by_cycle.pdf"
#    File quality_by_cycle_metrics = "${output_bam_prefix}.quality_by_cycle_metrics"
#    File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
#    File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
#  }
}


# Mark duplicate reads to avoid counting non-independent observations
# Aggregate aligned+merged flowcell BAM files and mark duplicates
# We take advantage of the tool's ability to take multiple BAM inputs and write out a single output
# to avoid having to spend time just merging BAM files.
function MarkDuplicates() {
    if [ $# -lt 1 ] ; then 
        echo "Err: ${FUNCNAME[0]} [output_base_name] input_bam ..." ; exit 1 
    else echo "${FUNCNAME[0]}" $* ; fi
    ###
    ###	CAUTION CAUTION CAUTION CAUTION CAUTION
    ### INCONSISTENT USAGE
    ###
    ### If there is only one argument it is a bam file to process, output
    ### will be derived from the bam file name
    ###
    ### If there is more than one, the first argument is the output base name
    ###
    ###
    if [ $# -gt 1 ] ; then
        # first argument is the base name of the output
        # following arguments are the bam file(s) to process
        ensemble_name=$1
        shift
    else
        ensemble_name=""
    fi
    #  Array[File] input_bams
    #  String output_bam_basename
    #  String metrics_filename
    #  external Int compression_level
    #  String? read_name_regex
    input_bams=( $@ )
    read_name_regex=${read_name_regex:-""}
    
    if [ ${#input_bams[@]} -gt 1 -o "$ensemble_name" != "" ] ; then
        # use a summary name
        #ensemble_name=$study_name.$ref_name	# defined at the top
        output_bam_basename=${ensemble_name}
    else
        # use a file specific name
        input_bam=$input_bams
        # remove trailing .bam
        output_bam_basename=${input_bam%.bam}
    fi
    output_bam=${output_bam_basename}.duplicates_marked.bam
    metrics_filename=${output_bam_basename}.duplicate_metrics

    # prepare argument list of input files
    input_list=''
    for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done
    
    
    if [ ! -s ${output_bam} ] ; then
      # The program default for READ_NAME_REGEX is appropriate in nearly every case.
      # Sometimes we wish to supply "null" in order to turn off optical duplicate detection
      # This can be desirable if you don't mind the estimated library size being wrong and optical duplicate detection is taking >7 days and failing

      # Task is assuming query-sorted input so that the Secondary and Supplementary reads get marked correctly
      # This works because the output of BWA is query-grouped and therefore, so is the output of MergeBamAlignment.
      # While query-grouped isn't actually query-sorted, it's good enough for MarkDuplicates with ASSUME_SORT_ORDER="queryname"
      ${picard_exec} \
            --java-options "-Dsamjdk.compression_level=${compression_level} \
              ${MarkDuplicates_java_opt}" \
            MarkDuplicates \
            ${input_list} \
            --OUTPUT ${output_bam} \
            --METRICS_FILE ${metrics_filename} \
            --VALIDATION_STRINGENCY SILENT \
            --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 \
            --ASSUME_SORT_ORDER "queryname" \
            --CREATE_MD5_FILE true \
            --CLEAR_DT "false" \
            --ADD_PG_TAG_TO_READS "false" \
            --READ_NAME_REGEX \
            #--REMOVE_DUPLICATES "true"

    fi
    #  output {
    #    File output_bam = "${output_bam_basename}.bam"
    #    File duplicate_metrics = "${metrics_filename}"
    #  }
}


# Sort BAM file by coordinate order and fix tag values for NM and UQ
function SortSam() {
  if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_bam
  #String output_bam_basename
  input_bam=$1
  output_bam_basename=${input_bam/%.bam/.sorted}

  if [ ! -s ${output_bam_basename}.bam ] then
    ${picard_exec} \
      --java-options "-Dsamjdk.compression_level=${compression_level} \
        -Xms4000m" \
      SortSam \
      --INPUT ${input_bam} \
      --OUTPUT ${output_bam_basename}.bam \
      --SORT_ORDER "coordinate" \
      --CREATE_INDEX "true" \
      --CREATE_MD5_FILE "true" \
      --MAX_RECORDS_IN_RAM 300000

  fi
  #output {
  #  File output_bam = "${output_bam_basename}.bam"
  #  File output_bam_index = "${output_bam_basename}.bai"
  #  File output_bam_md5 = "${output_bam_basename}.bam.md5"
  #}
}

# Collect alignment summary and GC bias quality metrics
function CollectReadgroupBamQualityMetrics () 
{
  if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_bam
  #File input_bam_index
  #String output_bam_prefix
  #Global File ref_dict
  #Global File ref_fasta
  #Global File ref_fasta_index
  input_bam=$1
  output_bam_prefix=${input_bam/%.bam/.quality}
  
  if [ ! -e "${output_bam_prefix}.alignment_summary_metrics" ] ; then
    ${picard_exec} \
        --java-options "-Xms5000m" \
        CollectMultipleMetrics \
        --INPUT ${input_bam} \
        --REFERENCE_SEQUENCE ${ref_fasta} \
        --OUTPUT ${output_bam_prefix} \
        --ASSUME_SORTED "true" \
        --PROGRAM "null" \
        --PROGRAM "CollectAlignmentSummaryMetrics" \
        --PROGRAM "CollectGcBiasMetrics" \
        --METRIC_ACCUMULATION_LEVEL "null" \
        --METRIC_ACCUMULATION_LEVEL "READ_GROUP"
  fi
  #output {
  #  File alignment_summary_metrics = "${output_bam_prefix}.alignment_summary_metrics"
  #  File gc_bias_detail_metrics = "${output_bam_prefix}.gc_bias.detail_metrics"
  #  File gc_bias_pdf = "${output_bam_prefix}.gc_bias.pdf"
  #  File gc_bias_summary_metrics = "${output_bam_prefix}.gc_bias.summary_metrics"
  #}
}

# Collect quality metrics from the aggregated bam
function CollectAggregationMetrics ()
{
  if [ $# -ne 1 ] ; then echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_bam
  #File input_bam_index
  #String output_bam_prefix
  #Global File ref_dict
  #Global File ref_fasta
  #Global File ref_fasta_index
  input_bam=$1
  output_bam_prefix=${input_bam/%.bam/.aggregate_quality}

  if [ ! -e "${output_bam_prefix}.alignment_summary_metrics" ] ; then
    ${picard_exec} \
      --java-options "-Xms5000m" \
      CollectMultipleMetrics \
      --INPUT ${input_bam} \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --OUTPUT ${output_bam_prefix} \
      --ASSUME_SORTED true \
      --PROGRAM "null" \
      --PROGRAM "CollectAlignmentSummaryMetrics" \
      --PROGRAM "CollectInsertSizeMetrics" \
      --PROGRAM "CollectSequencingArtifactMetrics" \
      --PROGRAM "CollectGcBiasMetrics" \
      --PROGRAM "QualityScoreDistribution" \
      --METRIC_ACCUMULATION_LEVEL "null" \
      --METRIC_ACCUMULATION_LEVEL "SAMPLE" \
      --METRIC_ACCUMULATION_LEVEL "LIBRARY"

    touch ${output_bam_prefix}.insert_size_metrics
    touch ${output_bam_prefix}.insert_size_histogram.pdf
  fi
  
  #output {
  #  File alignment_summary_metrics = "${output_bam_prefix}.alignment_summary_metrics"
  #  File bait_bias_detail_metrics = "${output_bam_prefix}.bait_bias_detail_metrics"
  #  File bait_bias_summary_metrics = "${output_bam_prefix}.bait_bias_summary_metrics"
  #  File gc_bias_detail_metrics = "${output_bam_prefix}.gc_bias.detail_metrics"
  #  File gc_bias_pdf = "${output_bam_prefix}.gc_bias.pdf"
  #  File gc_bias_summary_metrics = "${output_bam_prefix}.gc_bias.summary_metrics"
  #  File insert_size_histogram_pdf = "${output_bam_prefix}.insert_size_histogram.pdf"
  #  File insert_size_metrics = "${output_bam_prefix}.insert_size_metrics"
  #  File pre_adapter_detail_metrics = "${output_bam_prefix}.pre_adapter_detail_metrics"
  #  File pre_adapter_summary_metrics = "${output_bam_prefix}.pre_adapter_summary_metrics"
  #  File quality_distribution_pdf = "${output_bam_prefix}.quality_distribution.pdf"
  #  File quality_distribution_metrics = "${output_bam_prefix}.quality_distribution_metrics"
  #}

}

# Check that the fingerprints of separate readgroups all match
function CrossCheckFingerprints() {
    if [ $# -lt 1 ] ; then 
        echo "Err: ${FUNCNAME[0]} [output_base_name] input_bam ..." ; exit 1 
    else echo "${FUNCNAME[0]}" $* ; fi
    ###
    ###	CAUTION CAUTION CAUTION CAUTION CAUTION
    ### INCONSISTENT USAGE
    ###
    ### If there is only one argument it is a bam file to process, output
    ### will be derived from the bam file name
    ###
    ### If there is more than one, the first argument is the output base name
    ###
    ###
    if [ $# -gt 1 ] ; then
        # first argument is the base name of the output
        # following arguments are the bam file(s) to process
        ensemble_name=$1
        shift
    else
        ensemble_name=""
    fi
    #Array[File] input_bams
    #Array[File] input_bam_indexes
    #global File? haplotype_database_file
    input_bams=( $@ )
    
    if [ ${#input_bams[@]} -gt 1 -o "$ensemble_name" != "" ] ; then
        # use a summary name
        #ensemble_name=$study_name.$ref_name	# defined at the top
        output_basename=${ensemble_name}
    else
        # use a file specific name
        input_bam=$input_bams
        # remove trailing .bam
        output_basename=${input_bam%.bam}
    fi
    metrics_filename=${output_basename}.crosscheck
  
    # prepare argument list of input files
    input_list=''
    for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done

    if [ ! -e ${metrics_filename} ] ; then
      ${picard_exec} \
        --java-options "-Dsamjdk.buffer_size=131072 -XX:GCTimeLimit=50 \
          -XX:GCHeapFreeLimit=10 -Xms2000m" \
        CrosscheckReadGroupFingerprints \
        $input_list \
        --OUTPUT ${metrics_filename} \
        --HAPLOTYPE_MAP ${haplotype_database_file} \
        --EXPECT_ALL_READ_GROUPS_TO_MATCH "true" \
        --LOD_THRESHOLD -20.0
    fi
    #output {
    #  File metrics = "${metrics_filename}"
}

# Check that the fingerprint of the sample BAM matches the sample array
function CheckFingerprint() {
  if [ $# -ne 2 ] ; then echo "Err: ${FUNCNAME[0]} input_bam sample" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  
  #File input_bam
  #File input_bam_index
  #String output_basename
  #global File? haplotype_database_file
  #global File? genotypes
  #global File? genotypes_index
  #String sample
  input_bam=$1
  sample=$2
  output_basename=${input_bam%.bam}

  if [ ! -e "${output_basename}.fingerprinting_summary_metrics" ] ; then
    ${picard_exec} -Dsamjdk.buffer_size=131072 \
      -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms1024m  \
      CheckFingerprint \
      --INPUT ${input_bam} \
      --OUTPUT ${output_basename} \
      --GENOTYPES ${genotypes} \
      --HAPLOTYPE_MAP ${haplotype_database_file} \
      --SAMPLE_ALIAS "${sample}" \
      --IGNORE_READ_GROUPS "true"

  fi
  #output 
  #  File summary_metrics = "${output_basename}.fingerprinting_summary_metrics"
  #  File detail_metrics = "${output_basename}.fingerprinting_detail_metrics"
}

# Generate sets of intervals for scatter-gathering over chromosomes
function CreateSequenceGroupingTSV() {
  if [ $# -ne 0 ] ; then 
      echo "Err: ${FUNCNAME[0]}" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  #global File ref_dict

  # Use python to create the Sequencing Groupings used for BQSR and PrintReads Scatter.
  # It outputs to stdout where it is parsed into a wdl Array[Array[String]]
  # e.g. [["1"], ["2"], ["3", "4"], ["5"], ["6", "7", "8"]]
  if [ ! -e "sequence_grouping.txt" ] ; then
  (
    python <<CODE
    with open("${ref_dict}", "r") as ref_dict_file:
        sequence_tuple_list = []
        longest_sequence = 0
        for line in ref_dict_file:
            if line.startswith("@SQ"):
                line_split = line.split("\t")
                # (Sequence_Name, Sequence_Length)
                sequence_tuple_list.append((line_split[1].split("SN:")[1], int(line_split[2].split("LN:")[1])))
        longest_sequence = sorted(sequence_tuple_list, key=lambda x: x[1], reverse=True)[0][1]
    # We are adding this to the intervals because hg38 has contigs named with embedded colons and a bug in GATK strips off
    # the last element after a :, so we add this as a sacrificial element.
    hg38_protection_tag = ":1+"
    # initialize the tsv string with the first sequence
    tsv_string = sequence_tuple_list[0][0] + hg38_protection_tag
    temp_size = sequence_tuple_list[0][1]
    for sequence_tuple in sequence_tuple_list[1:]:
        if temp_size + sequence_tuple[1] <= longest_sequence:
            temp_size += sequence_tuple[1]
            tsv_string += "\t" + sequence_tuple[0] + hg38_protection_tag
        else:
            tsv_string += "\n" + sequence_tuple[0] + hg38_protection_tag
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
  )
  fi
  
  #output {
  #  Array[Array[String]] sequence_grouping = read_tsv("sequence_grouping.txt")
  #  Array[Array[String]] sequence_grouping_with_unmapped = read_tsv("sequence_grouping_with_unmapped.txt")

}

# Generate Base Quality Score Recalibration (BQSR) model
function BaseRecalibrator() {
    if [ $# -ne 2 ] ; then echo "Err: ${FUNCNAME[0]} input_bam intervals" ; exit 1 
    else echo "${FUNCNAME[0]}" $* ; fi

    #String input_bam
    #String recalibration_report_filename
    #Array[String] sequence_group_interval
    #global File dbSNP_vcf
    #global File dbSNP_vcf_index
    #global Array[File] known_indels_sites_VCFs
    #global Array[File] known_indels_sites_indices
    #global File ref_dict
    #global File ref_fasta
    #global File ref_fasta_index
    input_bam=$1
    shift
    recalibration_report_filename=${input_bam%.bam}.recal_data.csv
    sequence_group_interval=( $@ )
    
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

    $gatk_exec \
      --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -XX:+PrintFlagsFinal \
        -XX:+PrintGCTimeStamps -XX:+PrintGCDateStamps -XX:+PrintGCDetails \
        -Xloggc:gc_log.log -Xms4000m" \
      BaseRecalibrator \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${recalibration_report_filename} \
      -knownSites ${dbSNP_vcf} \
      $known_indels_arg \
      $sg_arg
  }

  #output {
  #  File recalibration_report = "${recalibration_report_filename}"
  #}
}



# Apply Base Quality Score Recalibration (BQSR) model
function ApplyBQSR() {
  if [ $# -ne 3 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam bqsr_report intervals" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #String input_bam
  #String output_bam_basename
  #File recalibration_report
  #Array[String] sequence_group_interval
  #global File ref_dict
  #global File ref_fasta
  #global File ref_fasta_index
  #global Int compression_level
  input_bam=$1
  shift
  output_bam_basename=${input_bam%.bam}.recalibrated
  recalibration_report=$2
  shift
  sequence_group_interval=( $@ )
  
  sg_arg=''
  for i in ${sequence_group_interval[@]} ; do
      sg_arg="${sg_arg} -L $i"
  done

  output_bam=$output_bam_basename.bam

  if [ ! -e $output_bam ] ; then
    $gatk_exec \
      --java-options "-XX:+PrintFlagsFinal -XX:+PrintGCTimeStamps \
        -XX:+PrintGCDateStamps -XX:+PrintGCDetails -Xloggc:gc_log.log \
        -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 \
        -Dsamjdk.compression_level=${compression_level} -Xms3000m" \
      ApplyBQSR \
      --createOutputBamMD5 \
      --add-output-sam-program-record \
      -R ${ref_fasta} \
      -I ${input_bam} \
      --useOriginalQualities \
      -O ${output_bam_basename}.bam \
      -bqsr ${recalibration_report} \
      -SQQ 10 -SQQ 20 -SQQ 30 \
      $sg_arg
  fi
  #output {
  #  File recalibrated_bam = "${output_bam_basename}.bam"
  #  File recalibrated_bam_checksum = "${output_bam_basename}.bam.md5"
  #}
}

# Combine multiple recalibration tables from scattered BaseRecalibrator runs
function GatherBqsrReports() {
  if [ $# -lt 3 ] ; then 
      echo "Err: ${FUNCNAME[0]} ensemble_name input_bam input_bam ..." ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

    ensemble_name=$1
    shift
    input_bqsr_reports=( $@ )
  

    inrep_arg=''
    for i in ${input_bqsr_reports[@]} ; do
        inrep_arg="${inrep_arg} -I $i"
    done
    output_report_filename=${ensemble_name}.all.recal_data.csv

    if [ ! -s ${output_report_filename} ] ; then
        ${gatk_exec} --java-options "${GatherBqsrReports_java_opt}" \
            GatherBQSRReports \
            $inrep_arg \
            -O ${output_report_filename}
    fi
  #output {
  #  File output_bqsr_report = "${output_report_filename}"
  #}
}

# Combine multiple recalibrated BAM files from scattered ApplyRecalibration runs
function GatherBamFiles() {
  if [ $# -lt 3 ] ; then 
      echo "Err: ${FUNCNAME[0]} ensemble_name input_bam input_bam ..." ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  
    ensemble_name=$1
    shift
    input_bams=( $@ )
  
    in_bams_arg=''
    for i in ${input_bams[@]} ; do  
        in_bams_arg="${in_bams_arg} --INPUT $i"
    done
    output_bam=$ensemble_name.all.recal.bam

    if [ ! -s $output_bam ] ; then
        ${gatk_exec} --java-options "-Dsamjdk.compression_level=${compression_level} ${GatherBamFiles_java_opt}" \
             GatherBamFiles \
             ${in_bams_arg} \
             --OUTPUT ${output_bam} \
             --CREATE_INDEX true \
             --CREATE_MD5_FILE true
    fi
  #output {
  #  File output_bam = "${output_bam_basename}.bam"
  #  File output_bam_index = "${output_bam_basename}.bai"
  #  File output_bam_md5 = "${output_bam_basename}.bam.md5"
  #}
}

function CheckPreValidation() {
  if [ $# -lt 4 ] ; then 
      echo "Err: ${FUNCNAME[0]} dup_metrics chim_metrics maxdup maxchim" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File duplication_metrics
  #File chimerism_metrics
  #Float max_duplication_in_reasonable_sample
  #Float max_chimerism_in_reasonable_sample
  duplication_metrics=$1
  chimerism_metrics=$2
  max_duplication_in_reasonable_sample=$3
  max_chimerism_in_reasonable_sample=$4

  if [ ! -e "duplication_value.txt" -o ! -e "chimerism_value.txt" ] ; then
  (
    set -o pipefail
    set -e

    grep -A 1 PERCENT_DUPLICATION ${duplication_metrics} > duplication.csv
    grep -A 3 PCT_CHIMERAS ${chimerism_metrics} | grep -v OF_PAIR > chimerism.csv

    python <<CODE

    import csv
    with open('duplication.csv') as dupfile:
      reader = csv.DictReader(dupfile, delimiter='\t')
      for row in reader:
        with open("duplication_value.txt","w") as file:
          file.write(row['PERCENT_DUPLICATION'])
          file.close()

    with open('chimerism.csv') as chimfile:
      reader = csv.DictReader(chimfile, delimiter='\t')
      for row in reader:
        with open("chimerism_value.txt","w") as file:
          file.write(row['PCT_CHIMERAS'])
          file.close()

    CODE

  )
  fi
  
  #output {
  #  Float duplication_rate = read_float("duplication_value.txt")
  #  Float chimerism_rate = read_float("chimerism_value.txt")
  #  Boolean is_outlier_data = duplication_rate > max_duplication_in_reasonable_sample || chimerism_rate > max_chimerism_in_reasonable_sample
  #}
}

function ValidateSamFile {
  if [ $# -lt 4 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_cram is_outlier max_out [ignore_str]" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
#  File input_bam
#  File? input_bam_index
#  String report_filename
#  File ref_dict
#  File ref_fasta
#  File ref_fasta_index
#  Int? max_output
#  Array[String]? ignore
#  Boolean? is_outlier_data
  input_bam=$1
  is_outlier_data=$2
  max_output=$3
  ignore=$4
  report_filename=${input_bam}.validation_report
  
  if [ ! -e $report_filename ] ; then
    ${picard_exec} \
      --java-options " -Xms6000m" \
      ValidateSamFile \
      --INPUT ${input_bam} \
      --OUTPUT ${report_filename} \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --MAX_OUTPUT ${max_output} \
      --IGNORE ${ignore:-"null"} \
      --MODE VERBOSE \
      --SKIP_MATE_VALIDATION ${is_outlier_data:-"false"} \
      --IS_BISULFITE_SEQUENCED false
  }
  #output {
  #  File report = "${report_filename}"
  #}
}


# IMPORTANT
# Note that these tasks will break if the read lengths in the bam are greater 
# than 250.
function CollectWgsMetrics() {
  if [ $# -ne 3 -o $# -ne 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam interval_list [read_size]" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

#  File input_bam
#  File input_bam_index
#  String metrics_filename
#  File wgs_coverage_interval_list
#  global File ref_fasta
#  global File ref_fasta_index
#  Int? read_length
  input_bam=$1
  metrics_filename=${input_bam%.bam}.wgs_metrics
  wgs_coverage_interval_list=$2
  read_length=${3:-250}

  if [ ! -e $metrics_filename ] ; then
    ${picard_exec} \
      --java-options " -Xms2000m" \
      CollectWgsMetrics \
      --INPUT ${input_bam} \
      --VALIDATION_STRINGENCY "SILENT" \
      --REFERENCE_SEQUENCE=${ref_fasta} \
      --INCLUDE_BQ_HISTOGRAM "true" \
      --INTERVALS ${wgs_coverage_interval_list} \
      --OUTPUT ${metrics_filename} \
      --USE_FAST_ALGORITHM true \
      --READ_LENGTH ${read_length:-250}
  fi
  #output {
  #  File metrics = "${metrics_filename}"
  #}
}

# Collect raw WGS metrics (commonly used QC thresholds)
function CollectRawWgsMetrics {
  if [ $# -ne 3 -o $# -ne 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam interval_list [read_size]" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

#  File input_bam
#  File input_bam_index
#  String metrics_filename
#  File wgs_coverage_interval_list
#  global File ref_fasta
#  global File ref_fasta_index
#  Int? read_length
  input_bam=$1
  metrics_filename=${input_bam%.bam}.raw_wgs_metrics
  wgs_coverage_interval_list=$2
  read_length=${3:-250}

  if [ ! -e $metrics_filename ] ; then
    ${picard_exec} \
      --java-options "-Xms2000m" \
      CollectRawWgsMetrics \
      --INPUT ${input_bam} \
      --VALIDATION_STRINGENCY SILENT \
      --REFERENCE_SEQUENCE ${ref_fasta} \
      --INCLUDE_BQ_HISTOGRAM true \
      --INTERVALS ${wgs_coverage_interval_list} \
      --OUTPUT ${metrics_filename} \
      --USE_FAST_ALGORITHM true \
      --READ_LENGTH ${read_length:-250}
  fi
  #output {
  #  File metrics = "${metrics_filename}"
  #}
}

# Generate a checksum per readgroup
function CalculateReadGroupChecksum() {
  if [ $# -ne 1 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_bam
  #File input_bam_index
  #String read_group_md5_filename
  
  read_group_md5_filename=${input_bam}.read_group_md5
  
  if [ ! -e "read_group_md5_filename" ] ; then
    ${picard_exec} \
      --java-options "-Xms1000m" \
      CalculateReadGroupChecksum \
      --INPUT ${input_bam} \
      --OUTPUT ${read_group_md5_filename}
  fi
  #output {
  #  File md5_file = "${read_group_md5_filename}"
  #}
}

# Notes on the contamination estimate:
# The contamination value is read from the FREEMIX field of the selfSM file output by verifyBamId
#
# In Zamboni production, this value is stored directly in METRICS.AGGREGATION_CONTAM
#
# Contamination is also stored in GVCF_CALLING and thereby passed to HAPLOTYPE_CALLER
# But first, it is divided by an underestimation factor thusly:
#   float(FREEMIX) / ContaminationUnderestimationFactor
#     where the denominator is hardcoded in Zamboni:
#     val ContaminationUnderestimationFactor = 0.75f
#
# Here, I am handling this by returning both the original selfSM file for reporting, and the adjusted
# contamination estimate for use in variant calling
#
# NOTE: This function produces a floating number in stdout and a selfSM file.
#
function CheckContamination() {
  if [ $# -ne 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam underestimation_factor" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  
  #File input_bam
  #File input_bam_index
  #global File contamination_sites_ud
  #global File contamination_sites_bed
  #global File contamination_sites_mu
  #global File ref_fasta
  #global File ref_fasta_index
  #String output_prefix
  #Float contamination_underestimation_factor
  input_bam=$1
  contamination_underestimation_factor=$2
  
  output_prefix=${input_bam%.bam}.preBqsr

  if [ ! -e ${output_prefix}.selfSM ] ; then
  (
    set -e

    # creates a ${output_prefix}.selfSM file, a TSV file with 2 rows, 19 columns.
    # First row are the keys (e.g., SEQ_SM, RG, FREEMIX), second row are the associated values
    $VerifyBamID_exec \
    --Verbose \
    --NumPC 4 \
    --Output ${output_prefix} \
    --BamFile ${input_bam} \
    --Reference ${ref_fasta} \
    --UDPath ${contamination_sites_ud} \
    --MeanPath ${contamination_sites_mu} \
    --BedPath ${contamination_sites_bed} \
    1>/dev/null

    # used to read from the selfSM file and calculate contamination, which gets printed out
    python3 <<CODE
    import csv
    import sys
    with open('${output_prefix}.selfSM') as selfSM:
      reader = csv.DictReader(selfSM, delimiter='\t')
      i = 0
      for row in reader:
        if float(row["FREELK0"])==0 and float(row["FREELK1"])==0:
          # a zero value for the likelihoods implies no data. This usually indicates a problem rather than a real event.
          # if the bam isn't really empty, this is probably due to the use of a incompatible reference build between
          # vcf and bam.
          sys.stderr.write("Found zero likelihoods. Bam is either very-very shallow, or aligned to the wrong reference (relative to the vcf).")
          sys.exit(1)
        print(float(row["FREEMIX"])/${contamination_underestimation_factor})
        i = i + 1
        # there should be exactly one row, and if this isn't the case the format of the output is unexpectedly different
        # and the results are not reliable.
        if i != 1:
          sys.stderr.write("Found %d rows in .selfSM file. Was expecting exactly 1. This is an error"%(i))
          sys.exit(2)
    CODE
  )
  fi
  #output {
  #  File selfSM = "${output_prefix}.selfSM"
  #  Float contamination = read_float(stdout())
}


# This task calls picard's IntervalListTools to scatter the input interval list into scatter_count sub interval lists
# Note that the number of sub interval lists may not be exactly equal to scatter_count.  There may be slightly more or less.
# Thus we have the block of python to count the number of generated sub interval lists.
function ScatterIntervalList() {
  if [ $# -ne 3 ] ; then 
      echo "Err: ${FUNCNAME[0]} interval_list scatter_count mult_break" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  #File interval_list
  #Int scatter_count
  #Int break_bands_at_multiples_of
  interval_list=$1
  scatter_count=$2
  break_bands_at_multiples_of=$3

  if [ ! -d 'out' ] ; then
  (
    set -e
    mkdir out
    ${picard_exec}
      --java-options " -Xms1g" \
      IntervalListTools \
      --SCATTER_COUNT ${scatter_count} \
      --SUBDIVISION_MODE BALANCING_WITHOUT_INTERVAL_SUBDIVISION_WITH_OVERFLOW \
      --UNIQUE true \
      --SORT true \
      --BREAK_BANDS_AT_MULTIPLES_OF ${break_bands_at_multiples_of} \
      --INPUT ${interval_list} \
      --OUTPUT out

    python3 <<CODE
    import glob, os
    # Works around a JES limitation where multiples files with the same name overwrite each other when globbed
    intervals = sorted(glob.glob("out/*/*.interval_list"))
    for i, interval in enumerate(intervals):
      (directory, filename) = os.path.split(interval)
      newName = os.path.join(directory, str(i + 1) + filename)
      os.rename(interval, newName)
    print(len(intervals))
    CODE
  )
  fi
  
  #output {
  #  Array[File] out = glob("out/*/*.interval_list")
  #  Int interval_count = read_int(stdout())
  #}
}

# Call variants on a single sample with HaplotypeCaller to produce a GVCF
function HaplotypeCaller() {
  if [ $# -ne 3 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam interval_list [contamination]" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #String input_bam
  #File interval_list
  #String gvcf_basename
  #global File ref_dict
  #global File ref_fasta
  #global File ref_fasta_index
  #Float? contamination
  input_bam=$1
  interval_list=$2
  contamination=$3
  gvcf_basename=${input_bam%.bam}

  # We use interval_padding 500 below to make sure that the HaplotypeCaller has context on both sides around
  # the interval because the assembly uses them.
  #
  # Using PrintReads is a temporary solution until we update HaploypeCaller to use GATK4. Once that is done,
  # HaplotypeCaller can stream the required intervals directly from the cloud.
  #
  # We are using GATK4, so the first step should not be needed
  if [ ! -e ${gvcf_basename}.vcf.gz ] ; then
    ${gatk_exec} \
      --java-options "-Xms2g" \
      PrintReads \
      -I ${input_bam} \
      --interval_padding 500 \
      -L ${interval_list} \
      -O ${gvcf_basename}.reads.bam \
    && \
    ${gatk_exec} \
      --java-options "-XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10 -Xms8000m" \
      HaplotypeCaller \
      -R ${ref_fasta} \
      -O ${gvcf_basename}.vcf.gz \
      -I  ${gvcf_basename}.reads.bam\
      -L ${interval_list} \
      -ERC GVCF \
      --max_alternate_alleles 3 \
      -contamination ${contamination:-0} \
      --read_filter OverclippedRead \
      #-bamout ${gvcf_basename}.reads.HC.bam \
      #-variant_index_parameter 128000 \
      #-variant_index_type LINEAR \
  fi
  #output {
  #  File output_gvcf = "${gvcf_basename}.vcf.gz"
  #  File output_gvcf_index = "${gvcf_basename}.vcf.gz.tbi"
  #}
}

# Combine multiple VCFs or GVCFs from scattered HaplotypeCaller runs
function MergeVCFs() {
  if [ $# -lt 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} [output_vcf_name] input_bam ..." ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #Array[File] input_vcfs
  #Array[File] input_vcfs_indexes
  #String output_vcf_name
  output_vcf_name=$1
  shift
  input_vcfs=( $@ )

  # Using MergeVcfs instead of GatherVcfs so we can create indices
  # See https://github.com/broadinstitute/picard/issues/789 for relevant GatherVcfs ticket
  # prepare argument list of input files
  input_list=''
  for i in ${input_bams[@]} ; do input_list="$input_list --INPUT $i" ; done

  if [ ! -e $output_vcf_name ] ; then
    ${picard_exec} \
      --java-options "-Xms2000m" \
      MergeVcfs \
      ${input_list} \
      --OUTPUT ${output_vcf_name}
  }
  #output {
  #  File output_vcf = "${output_vcf_name}"
  #  File output_vcf_index = "${output_vcf_name}.tbi"
  #}
}


# Validate a GVCF with -gvcf specific validation
function ValidateGVCF() {
  if [ $# -lt 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_gvcf interval_list" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi
  #File input_vcf
  #File input_vcf_index
  #global File ref_fasta
  #global File ref_fasta_index
  #global File ref_dict
  #global File dbSNP_vcf
  #global File dbSNP_vcf_index
  #global File wgs_calling_interval_list
  input_vcf=$1
  wgs_calling_interval_list=$2


    # NOTE: this may fail validation on exome data
    # which is why we use --warnOnErrors
    ${gatk_exec} \
       --java-options "-Xms3000m" \
      ValidateVariants \
      -V ${input_vcf} \
      -R ${ref_fasta} \
      -L ${wgs_calling_interval_list} \
      -gvcf \
      --validationTypeToExclude ALLELES \
      --dbsnp ${dbSNP_vcf} \
      --warnOnErrors
  
}

# Collect variant calling metrics from GVCF output
function CollectGvcfCallingMetrics() {
  if [ $# -lt 2 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_gvcf interval_list" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_vcf
  #File input_vcf_index
  #String metrics_basename
  #global File dbSNP_vcf
  #global File dbSNP_vcf_index
  #global File ref_dict
  #global File wgs_evaluation_interval_list
  input_vcf=$1
  wgs_evaluation_interval_list=$interval_list
  metrics_basename=${input_vcf%.vcf*}

  if [ ! -s ${output_basename}.variant_calling_summary_metrics \
    -o ! -s ${output_basename}.variant_calling_detail_metrics ] ; then
      ${picard_exec} \
        --java-options "-Xms2000m" \
        CollectVariantCallingMetrics \
        --INPUT ${input_vcf} \
        --OUTPUT ${metrics_basename} \
        --DBSNP ${dbSNP_vcf} \
        --SEQUENCE_DICTIONARY ${ref_dict} \
        --TARGET_INTERVALS ${wgs_evaluation_interval_list} \
        --GVCF_INPUT true
  fi
  #output {
  #  File summary_metrics = "${metrics_basename}.variant_calling_summary_metrics"
  #  File detail_metrics = "${metrics_basename}.variant_calling_detail_metrics"
  #}
}

# Convert BAM file to CRAM format
# Note that reading CRAMs directly with Picard is not yet supported
function ConvertToCram() {
  if [ $# -ne 1 ] ; then 
      echo "Err: ${FUNCNAME[0]} input_bam" ; exit 1 
  else echo "${FUNCNAME[0]}" $* ; fi

  #File input_bam
  #File ref_fasta
  #File ref_fasta_index
  #String output_basename
  input_bam=$1
  output_basename=${input_bam%.bam}

  if [ ! -e $output_basename.cram ] ; then
  (
    set -e
    set -o pipefail

    samtools view -C -T ${ref_fasta} ${input_bam} | \
    tee ${output_basename}.cram | \
    md5sum | awk '{print $1}' > ${output_basename}.cram.md5

    # Create REF_CACHE. Used when indexing a CRAM
    # this is a SAMtools script
    seq_cache_populate.pl -root ./ref/cache ${ref_fasta}
    export REF_PATH=:
    export REF_CACHE=./ref/cache/%2s/%2s/%s

    samtools index ${output_basename}.cram
  )
  
  #output {
  #  File output_cram = "${output_basename}.cram"
  #  File output_cram_index = "${output_basename}.cram.crai"
  #  File output_cram_md5 = "${output_basename}.cram.md5"
  #}
}

