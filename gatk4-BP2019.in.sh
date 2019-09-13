## SAMPLE NAME AND UNMAPPED BAMS
ref_name="hg19"		# genome
sample_name="PR1"	# sample
project_directory="${HOME}/work/lorena/proradium1"
flowcell_unmapped_bams_list="${project_directory}/uBam.list"
unmapped_bam_suffix=".bam"

## DATA DIRECTORY
data_dir="${HOME}/work/lorena/data/"
gatk_bundle_dir="${data_dir}/gatk-hg19-bundle"


## "PATHS",
bwa_path='/usr/bin/'
gatk_exec="${HOME}/contrib/gatk4/gatk"	# gatk executable
gotc_path="${HOME}/contrib/gatk4/"	# directory where programs reside
					# I suppose it refers to the directory
                                        # where GATK team have all their stuff
                                        # when working On The Cloud
