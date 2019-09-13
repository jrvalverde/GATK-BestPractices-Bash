# Assemble reads
#
# We assume that reads are mated pairs provided as two *_R[12]_001.fastq files
#
# We take care not to repeat any work, hence, re-running is safe.
#
# The graphics routines may fail on servers where no X support is available.
# If so, in order to re-run them, we need to remove the generated last, 
# mummer3 and nucmer directories and re-run in a different computer with
# appropriate support.
#
# If we want to genome-align against different reference genomes we can, either
# repeat the full process or only align by renaming reference.fna to something
# else (e.g. reference.fna.orig), renaming the last/mummer3/nucmer directories
# to allow re-running these programs against the new reference and 
# installing/symlinking a new reference. This is confusing, I know, and I
# should certainly separate assembly and comparison into separate programs
# one of these days, eventually.
#
# (C) José R. Valverde, CNB/CSIC. 2015-2016
#
set +x

# Default values
SEQUENCE_IS_454="no"		# -r --roche454
SEQUENCE_IS_ILLUMINA="no"	# -i --illumina

BWA="no"			# -b --BWA
BOWTIE2="no"			# -B --Bowtie2
LASTASS="no"			# -L --lastass
MAQ="no"			# -M --MAQ
A5="no"				# -A --A5
MUMMER="no"			# -m --mummer3
LAST="no"			# -l --last
NUCMER="no"			# -n --nucmer
VERBOSE=0                       # -v --verbose
SEQDIR="./orig"			# -S --sequences seqdir
R1="R1.fastq"			# -1 --R1
R2="R2.fastq"			# -2 --R2
REFFNA="./refs/reference.fna"	# -R --reference ref.fna

BOWTIE_DIR=~/contrib/bowtie2-2.2.4
export PATH=$BOWTIE_DIR:$PATH


function abspath() {
    # generate absolute path from relative path
    # $1     : relative filename
    # return : absolute path
    if [ -d "$1" ]; then
        # dir
        (cd "$1"; pwd)
    elif [ -f "$1" ]; then
        # file
        if [[ $1 == */* ]]; then
            echo "$(cd "${1%/*}"; pwd)/${1##*/}"
        else
            echo "$(pwd)/$1"
        fi
    else
        #echo "/dev/null"
        echo " ---FILE-DOES-NOT-EXIST---"
    fi
}

function relpath() {
    # both $1 and $2 are absolute paths beginning with /
    # returns relative path to $2($target) from $1($source)
    if [ $# -eq 2 ] ; then
       source="$1"
       target="$2"
    elif [ $# -eq 1 ] ; then
       source=`pwd`
       target="$1"
    else
       echo "---relpath error, usage: relpath [source] target---"
       return
    fi
    
    # example use (ensure a symlink contains the relative path to 
    # a file from our current position instead of its absolute path):
    #    f=`abspath $file`
    #    ln -sf `relpath $f` target.lnk
     
    if [ ! -e "$source" -o ! -e "$target" ] ; then
        echo "---NO-SUCH-FILE---"
        return
    fi

    # example:
    #    f=`abspath $file`
    #    ln -s `relpath `pwd` $f` f.lnk

    common_part=$source # for now
    result="" # for now

    while [[ "${target#$common_part}" == "${target}" ]]; do
        # no match, means that candidate common part is not correct
        # go up one level (reduce common part)
        common_part="$(dirname $common_part)"
        # and record that we went back, with correct / handling
        if [[ -z $result ]]; then
            result=".."
        else
            result="../$result"
        fi
    done

    if [[ $common_part == "/" ]]; then
        # special case for root (no common path)
        result="$result/"
    fi

    # since we now have identified the common part,
    # compute the non-common part
    forward_part="${target#$common_part}"

    # and now stick all parts together
    if [[ -n $result ]] && [[ -n $forward_part ]]; then
        result="$result$forward_part"
    elif [[ -n $forward_part ]]; then
        # extra slash removal
        result="${forward_part:1}"
    fi

    echo $result
}

function posix_relpath () {
    [ $# -ge 1 ] && [ $# -le 2 ] || return 1
    current="${2:+"$1"}"
    target="${2:-"$1"}"
    [ "$target" != . ] || target=/
    target="/${target##/}"
    [ "$current" != . ] || current=/
    current="${current:="/"}"
    current="/${current##/}"
    appendix="${target##/}"
    relative=''
    while appendix="${target#"$current"/}"
        [ "$current" != '/' ] && [ "$appendix" = "$target" ]; do
        if [ "$current" = "$appendix" ]; then
            relative="${relative:-.}"
            echo "${relative#/}"
            return 0
        fi
        current="${current%/*}"
        relative="$relative${relative:+/}.."
    done
    relative="$relative${relative:+${appendix:+/}}${appendix#/}"
    echo "$relative"
}

function fastasize() {
    awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' \
        "$1"
}

function fastasizes() {
    # adds total summaries at the end
    awk '/^>/ {if (seqlen){print seqlen};print;seqtotal+=seqlen;seqlen=0;seq+=1;next;}
        {seqlen=seqlen+length($0)}
        END{print seqlen;print seq" sequences, total length " seqtotal+seqlen}' \
        "$1"
}

function fq2f() {
    # Convert fastq file to fasta format
    # we'll simply strip the quality codes, without using them
    # output goes to stdout!!
    #set +x
    fq=`abspath "$1"`
    #echo "#fq2f: Converting $fq to fasta format"
    if [ ! -e "$fq" ] ; then
        echo "fq2f: error, nothing to do"
        return
    fi
    count=0
    deleted=0
    n=0
    # we have a problem when reading lines: by default read interprets 
    # backslashes, therefore, a line ending in \ will include the next 
    # line as well
    # We must use read -r to avoid this problem.
    cat "$fq" | \
    while read -r line ; do
        #n=$((n + 1))
        if [ "$line" = "+" ] ; then
            count=$((count - 1))	# discount the title line
            #echo "XXX $count $line"
            #while [ $count -gt 0 ] ; do
            while [ $deleted -lt $count ] ; do
                deleted=$((deleted + 1))
                read -r line
                #n=$((n + 1))
                #echo ">>>" $count $deleted $n $line
            done
            count=0	# reset counter
            deleted=0
        else
            #echo $count: $line
            echo $line
            count=$((count + 1))
         fi
    done | \
    tr '@' '>' 
    #set -x
}

function bwa_reference_assemble() {
    # assemble usinf a reference and R1 + R2 reads
    ref=`abspath "$1"`
    R1=`abspath "$2"`
    R2=`abspath "$3"`

    if [ ! -d bwa ] ; then
        echo "Creating bwa"
        mkdir bwa
    fi
    echo "Entering bwa"
    cd bwa
    if [ ! -e reference.fna ] ; then
        if [ -e "$ref" ] ; then
	    # this is safe because we know reference.fna does not exist
	    #	but it could be a dangling symlink, so we force removal
            #   of any existing symlink before creating the new one.
	    ln -sf `relpath "$ref"` reference.fna
            #ln -sf $ref reference.fna
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a reference sequence for BWA"
	    return
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
	    ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run BWA"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
	    ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run BWA"
            return
        fi
    fi

    # We require two mate-pair reads files in fastq format
    # and a reference file in fasta format in the current
    # directory
    #
    #	Further, they must be called reference.fna, R1.fastq and R2.fastq
    #
    # Check we have everything we need
    if [ ! -e reference.fna -o ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need a reference and two mate pairs files to run BWA"
        return
    fi
    
    # Run BWA
    # 1. Create reference database
    if [ ! -e reference.fna.bwt ] ; then
	echo "Making reference BWA database"
        bwa index -a is reference.fna
    fi
    # 2. Align
    if [ ! -e s1.bwa.sai ] ; then
	echo "Aligning R1"
        bwa aln -t 8 reference.fna R1.fastq > s1.bwa.sai
    fi
    if [ ! -e s2.bwa.sai ] ; then
	echo "Aligning R2"
        bwa aln -t 8 reference.fna R2.fastq > s2.bwa.sai
        # we may also try this
        #bwa mem -t 8 reference.fna R2.fastq > s2.bwa.sai
    fi
    if [ ! -e aligned.bwa.sam ] ; then
	echo "Aligning mate pairs"
        bwa sampe -P reference.fna \
    	    s1.bwa.sai s2.bwa.sai \
            R1.fastq R2.fastq > aligned.bwa.sam
    fi
    if [ -s aligned.bwa.cons.fastq -a ! -s aligned.bwa.cons.fna ] ; then
	echo "bwa_reference_assemble: making fasta from aligned.bwa.cons.fastq"
        fq2f aligned.bwa.cons.fastq > aligned.bwa.cons.fna
    else
        echo "bwa_reference_assemble: WARNING! Consensus fastq file size is zero"
    fi
    cd ..
}

function bowtie2_reference_assemble() {
    # assemble using a reference and R1, R2 and ($U)npaired reads
    ref=`abspath "$1"`
    R1=`abspath "$2"`
    R2=`abspath "$3"`
    U=`abspath "$4"`
    
    if [ ! -d bowtie2 ] ; then
        echo "Creating bowtie2"
        mkdir bowtie2
    fi
    echo "Entering bowtie"
    cd bowtie2
    if [ ! -e reference.fna ] ; then
        if [ -e "$ref" ] ; then
	    ln -sf `relpath "$ref"` reference.fna
            #ln -sf $ref reference.fna
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a reference sequence for Bowtie2"
	    return
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
	    ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run Bowtie2"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
	    ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run Bowtie2"
            return
        fi
    fi
    # U.fastq may be absent if there are no unpaired reads
    if [ -e "$U" ] ; then
        ln -sf `relpath "$U"` U.fastq
        #ln -s $U U.fastq
        # convert symlink from absolute to relative path using Mark Lord's
        #symlinks -c .
    fi
    # We require two mate-pair reads files in fastq format
    # and a reference file in fasta format in the current
    # directory
    #
    #	Further, they must be called reference.fna, R1.fastq and R2.fastq
    #
    if [ ! -e reference.fna -o ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need a reference and two mate pairs files to run Bowtie2"
        return
    fi

    # make reference database
    if [ ! -e reference.1.bt2 ] ; then
    	echo "Bowtie2: making reference database"
    	bowtie2-build reference.fna reference
    fi

    # align

    # check for availability of unaligned reads
    U=""
    if [ -e U.fastq ] ; then U="-U U.fastq" ; fi

    if [ ! -e aligned.bowtie.sam ] ; then
        echo "Bowtie2: aligning mate reads against reference"
        bowtie2 -p 8 --qc-filter -x reference \
	    -1 R1.fastq -2 R2.fastq $U\
            -S aligned.bowtie.sam
    fi
    cd ..
}


#
# map reads against a reference.
#	We need the reference in "reference.fna" and the reads (direct
#	and indirect) as "R1.fastq" and "R2.fastq"
#
last_reference_assemble() {
    # assemble usinf a reference and R1 + R2 reads
    ref=`abspath "$1"`
    R1=`abspath "$2"`
    R2=`abspath "$3"`

    if [ ! -d last ] ; then
        echo "Creating last"
        mkdir last
    fi
    echo "Entering last"
    cd last
    if [ ! -e reference.fna ] ; then
        if [ -e "$ref" ] ; then
	    # this is safe because we know reference.fna does not exist
	    #	but it could be a dangling symlink, so we force removal
            #   of any existing symlink before creating the new one.
	    ln -sf `relpath "$ref"` reference.fna
            #ln -sf $ref reference.fna
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a reference sequence for LAST"
	    return
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
	    ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run LAST"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
	    ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run LAST"
            return
        fi
    fi

    # We require two mate-pair reads files in fastq format
    # and a reference file in fasta format in the current
    # directory
    #
    #	Further, they must be called reference.fna, R1.fastq and R2.fastq
    #
    # Check we have everything we need
    if [ ! -e reference.fna -o ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need a reference and two mate pairs files to run LAST"
        return
    fi
    
    ##############
    # DO THE WORK
    ##############
    
    echo "Making reference database"
    lastdb ref reference.fna

    echo "Mapping individual reads"
    echo "R1"
    if [ ! -s R1.maf ] ; then
	lastal -Q1 -e120 -i1 ref R1.fastq > R1.maf
	last-split R1.maf > R1.split.maf
    fi
    echo "R2"
    if [ ! -s R2.maf ] ; then
	lastal -Q1 -e120 -i1 ref R2.fastq > R2.maf
	last-split R2.maf > R2.split.maf
    fi

    echo "Mapping paired reads"
    if [ ! -e paired.maf ] ; then
	last-pair-probs R1.maf R2.maf > paired.maf 2> last-pair-stats
    fi

    #last-pair-probs R1.split.maf R2.split.maf > split.mapped.maf 2> last-split-stats

    echo "Making SAM files and plots"
    for i in *.maf ; do
        echo "Converting $i to SAM and plotting"
	if [ ! -e `basename $i maf`sam ] ; then
            maf-convert sam $i > `basename $i maf`sam
	fi
        # This will likely fail because of too many reads
	if [ ! -e `basename $i maf`png ] ; then
            last-dotplot $i `basename $i maf`png
	fi
    done

    cd ..
}


function analyze_sam_file() {
    # IMPORTANT: FILE MUST BE IN CURRENT DIRECTORY!!!
    #	Usage: analyze_sam_file file.sam reference.fna
    
    # extract filename
    dir=`dirname "$1"`
    base=`basename "$1" .sam`
    cd "$dir"

    ref=`abspath "$2"`

    if [ ! -e "$1" -o ! -e "$2" ] ; then 
    	echo "analyze_sam_file: '$1' and '$2' must exist"
	return
    fi

    echo "Checking SAM file"
    head -n 1 $base.sam | grep '@SQ'
    if [ $? -eq 1 ] ; then
        echo "Sanitizing sam file (old file named $base.sam.orig)"
        # insert sequence information from the reference sequence
        mv $base.sam $base.sam.orig
        sq=`head -1 $ref | cut -d' ' -f1 | tr -d '>'`
	ln=`tail -n+2 $ref | tr -d '\n' | wc -c`
	echo "@SQ	SN:$sq	LN:$ln" | cat - $base.sam.orig > $base.sam
    fi

    echo "Analyzing SAM file"

    # Postprocess to BAM
    if [ ! -e ${base}.sorted.bam.bai ] ; then
	echo "Creating ${base}.sorted.bam.bai"
        samtools view -bS ${base}.sam -o ${base}.bam
        samtools sort ${base}.bam ${base}.sorted
        samtools index ${base}.sorted.bam
    fi

    # Do variant calling
    # variant calling (pre-samtools 1.x)
    #samtools faidx $ref
    #samtools pileup -f $ref ${base}.sorted.bam | \
    #	bcftools view -bvcg - > ${base}.bcf
    #bcftools view ${base}.bcf > ${base}.vcf
    if [ ! -e ${base}.variants.vcf ] ; then
	echo "Doing variant calling"
        samtools mpileup -d 10000 -g -f $ref ${base}.sorted.bam > \
 	    ${base}.variants.bcf
        bcftools view -c -v ${base}.variants.bcf > ${base}.variants.vcf
    fi

    vcfutils=/usr/share/samtools/vcfutils.pl
    # Do consensus
    # consensus (very old samtools) ***
    if [ ! -e ${base}.cons.fastq ] ; then
	echo "Calculating consensus"
	# old samtools
        #~/src/samtools/samtools-527/samtools pileup -cf $ref ${base}.sorted.bam > \
	#    ${base}.cons.pileup
        #perl /usr/share/samtools/samtools.pl pileup2fq ${base}.cons.pileup > \
	#    ${base}.cons.fastq
	
	# consensus (samtools 1.x) DOES NOW WORK
	samtools mpileup -uf $ref ${base}.sorted.bam | \
		bcftools view -gc - | $vcfutils vcf2fq > ${base}.cons.fastq
	#	bcftools view -gcv - | $vcfutils vcf2fq > ${base}.cons.fastq
	#	bcftools consensus -f $ref - | $vcfutils vcf2fq > ${base}.cons.fastq
	
	#samtools mpileup -d 10000 -v -f $ref ${base}.sorted.bam > E.vcf 
	#bcftools consensus -f $ref E.vcf -o E.c.vcf
	#cat E.c.vcf | vcfutils.pl vcf2fq > E.fq
    fi

    # Compute coverage
    if [ ! -e ${base}.cov.stats ] ; then
	echo "Computing coverage"
        samtools mpileup -d 10000 -ABR ${base}.sorted.bam > ${base}.coverage
        # requires processing to compute averages, mean, etc...
        cut -f4 ${base}.coverage | Rscript -e 'd<-scan("stdin", quiet=TRUE)' \
              -e 'cat(min(d), max(d), median(d), mean(d), "\n", sep="\t")' \
              -e 'summary(d)' > ${base}.cov.stats
    fi

    ###
    if [ ! -e ${base}.sorted.fastq ] ; then
	echo "Piling up"
        samtools mpileup -u -f $ref ${base}.sorted.bam > ${base}.sorted.bcf
        bcftools view -v -c -g ${base}.sorted.bcf > ${base}.sorted.vcf
        perl /usr/share/samtools/vcfutils.pl vcf2fq ${base}.sorted.vcf > ${base}.sorted.fastq
    fi
}

function test_reference_assembly() {
    ref=`abspath "$1"`
    contigs=`abspath "$2"`

}

function maq_reference_assemble() {
    ref=`abspath "$1"`
    R1=`abspath "$2"`
    R2=`abspath "$3"`

    if [ ! -d maq ] ; then
        echo "Creating maq"
        mkdir maq
    fi
    echo "Entering maq"
    cd maq
    if [ ! -e reference.fna ] ; then
        if [ -e "$ref" ] ; then
	    ln -sf `relpath "$ref"` reference.fna
            #ln -sf $ref reference.fna
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a reference sequence for MAQ"
	    return
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
	    ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run MAQ"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
            ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            # convert symlink from absolute to relative path using Mark Lord's
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run MAQ"
            return
        fi
    fi

    # We require two mate-pair reads files in fastq format
    # and a reference file in fasta format in the current
    # directory
    #
    #	Further, they must be called reference.fna, R1.fastq and R2.fastq
    #
    if [ ! -e reference.fna -o ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need a reference and two mate pairs files to run MAQ"
        return
    fi
    if [ ! -e easyrun.fq/cns.fq ] ; then
        echo "Reference-assembling using MAQ"
        maq.pl easyrun reference.fna R1.fastq R2.fastq > maq.out 2>&1
        # generates easyrun/cns.fq
        fq2f easyrun/cns.fq > easyrun/cns.fna
        mv easyrun easyrun.fq
    fi
    if [ ! -e easyrun.fa/cns.fq ] ; then
        # try again with fasta reads
        echo "MAQ: Converting R1 and R2 fastq files to fasta format for re-assembly"
        fq2f R1.fastq > R1.fa
        fq2f R2.fastq > R2.fa
        maq.pl easyrun reference.fna R1.fa R2.fa > maq.out 2>&1
        fq2f easyrun/cns.fq > easyrun/cns.fna
        mv easyrun easyrun.fa
    fi    

    cd ..
}

function a5_denovo_assemble() {
    R1=`abspath "$1"`
    R2=`abspath "$2"`

    if [ ! -d a5 ] ; then
        echo "Creating a5"
        mkdir a5
    fi
    echo "Entering a5"
    cd a5
    if [ ! -e R1.fastq ] ; then
        if [ -e "$R1" ] ; then
            ln -sf `relpath "$R1"` R1.fastq
            #ln -sf $R1 R1.fastq
            #symlinks -c .
        else
            echo "We need a filtered reads file (1) to run A5"
            return
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e "$R2" ] ; then
            ln -sf `relpath "$R2"` R2.fastq
            #ln -sf $R2 R2.fastq
            #symlinks -c .
        else
            echo "We need a filtered reads file (2) to run A5"
            return
        fi
    fi

    # We require two mate-pair reads files in the current directory
    # which must be called R1.fastq and R2.fastq
    #
    if [ ! -e R1.fastq -o ! -e R2.fastq ] ; then
        echo "We need two mate pairs files to run A5"
        return
    fi

    if [ ! -e genome.final.scaffolds.fastq ] ; then
	export PATH=~/contrib/a5/bin:$PATH
        a5_pipeline.pl R1.fastq R2.fastq genome
    fi
    cd ..
}

function mummer3_align() {
    if [ "$MUMMER" != "yes" ] ; then return ; fi
    ref=`abspath "$1"`
    qontigs=`abspath "$2"`

    #echo "mummer3 $ref $qontigs"    
    if [ ! -d mummer3 ] ; then
        echo "Creating mummer3"
        mkdir mummer3
    fi
    echo "Entering mummer3"
    cd mummer3
    if [ ! -e reference.fna ] ; then
        echo linking `relpath "$ref"` as reference.fna
	ln -sf `relpath "$ref"` reference.fna
        #ln -sf $ref reference.fna
        #symlinks -c .
    fi
    f="${qontigs##*/}"
    ext="${qontigs##*.}"
    name="${f%.*}"
    if [ ! -e contigs.fna ] ; then
        if [ "$ext" = "fastq" ] ; then
    	    # convert to fasta
            fq2f $qontigs > contigs.fna
        else
            ln -sf `relpath $qontigs` contigs.fna
            #ln -sf $qontigs contigs.fna
            #symlinks -c .
        fi
    fi
    out="$name"

    #set -x
    if [ ! -e $out.out ] ; then
        echo "Comparing $out contigs.fna against reference.fna"
        run-mummer3 reference.fna contigs.fna $out
	# plot
	# --oplor	color plot lines by percent similarity
	# -c		generate coverage plot
	# -f		only display best hit
	# -l            layout a .delta multiplot intelligibly
	# -R -Q
	# -t png
        mummerplot --color -c -t png -p ${out}-cov $out.out
	mummerplot            -t png -p ${out}-dp  $out.out
	
    fi
    
    if [ ! -e ${out}-concat.mums ] ; then
        echo "Concatenating contigs and aligning against refereece"
	echo ">concatenated_fragments" > ${out}-concat.fna
	grep -v "^>" contigs.fna >> ${out}-concat.fna
	mummer -mum -b -c reference.fna ${out}-concat.fna > ${out}-concat.mums
	mummerplot --color -c -t png -p ${out}-concat-cov ${out}-concat.mums
	mummerplot --color    -t png -p ${out}-concat-dp  ${out}-concat.mums
	mummer -mum -b -c ${out}-concat.fna reference.fna > ${out}-tacnoc.mums
	mummerplot --color -c -t png -p ${out}-tacnoc-cov ${out}-tacnoc.mums
	mummerplot --color    -t png -p ${out}-tacnoc-dp  ${out}-tacnoc.mums
    fi
    cd ..
    #set +x
}

function nucmer_align() {
    if [ "$NUCMER" != "yes" ] ; then return ; fi
    ref=`abspath $1`
    sqs=`abspath $2`
    
    if [ ! -e $1 -o ! -e $2 ] ; then
        echo "nucmer: $1 and $2 must exist"
        return
    fi
    #echo "nucmer $ref $qontigs"    
    if [ ! -d nucmer ] ; then
        echo "Creating nucmer"
        mkdir nucmer
    fi
    echo "Entering nucmer"
    cd nucmer
    if [ ! -e reference.fna ] ; then
        ln -sf `relpath $ref` reference.fna
        #ln -sf $ref reference.fna
        #symlinks -c .
    fi
    
    f="${sqs##*/}"
    ext="${sqs##*.}"
    name="${f%.*}"
    echo "nucmer_align: $sqs - $f - $name - $ext"
    if [ ! -e $name.fna ] ; then
        if [ "$ext" = "fastq" ] ; then
    	    # convert to fasta and use fasta file instead
	    echo "Converting $sqs to fasta format"
            fq2f $sqs > $name.fna
            sqs=$name.fna
        else
            ln -sf `relpath $sqs` $name.fna
            #ln -sf $sqs $name.fna
            #symlinks -c .
        fi
    fi
    out=$name

    #set -x
    #echo `awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' reference.fna`
    #reflen=`awk '/^>/ {if (seqlen){print seqlen}; print ;seqlen=0;next; } { seqlen = seqlen +length($0)}END{print seqlen}' \
    #        reference.fna | tail -1`
    reflen=`fastasize reference.fna | tail -1`
    echo "Referece" > SIZES
    fastasizes reference.fna >> SIZES
    echo "Probes" >> SIZES 
    fastasizes $name.fna >> SIZES

    if [ ! -s $out.delta ] ; then
        nucmer -maxmatch -p $out $ref $name.fna
	mummerplot --color -c -t png -p ${out}-cov $out.delta
	mummerplot --color    -t png -p ${out}-dp  $out.delta
    fi
    if [ ! -s $out.coords ] ; then
        show-coords -L 100 -I 50 -r -c -l $out.delta > $out.coords
	#mapview -n 1 -p ${out}-mapview $out.coords
	mapview -n 1 -f png -p ${out}-mapview $out.coords
	#mapview -n 1 -f pdf -p ${out}-mapview $out.coords
    fi
    
    if [ ! -s $out.snps ] ; then
        show-snps -C $out.delta > $out.snps
    fi
    
    if [ ! -s $out.aligns ] ; then
        show-aligns -r $out.delta $ref $sqs > $out.aligns
    fi
    
    if [ ! -s $out.tiling ] ; then
        # assume circular reference (-c)
        show-tiling -c -p $out.pseudo $out.delta > $out.tiling
	mummerplot --color -c -t png -p ${out}-tile-cov ${out}.tiling
    fi

    # compute coverage by summing-up inter-contig gaps
    if [ ! -s $out.tiling ] ; then
        echo "0" > $out.coverage
    else
        #cat $out.tiling | tail -n +2 | cut -d'	' -f3
	echo -n "coverage: "
        cat $out.tiling | tail -n +2 | cut -d'	' -f3 | \
	    awk "{s+=\$1} END {print (1-(s/$reflen))*100 }" 
        cat $out.tiling | tail -n +2 | cut -d'	' -f3 | \
	    awk "{s+=\$1} END {print (1-(s/$reflen))*100 }" > $out.coverage
    fi

    # Now, repeat for the concatenated contigs
    #
    if [ ! -e ${out}-concat.delta ] ; then
        echo "nucmer: Concatenating $sqs and comparing against reference"
        echo ">Concatenated_contigs" > ${out}-concat.fna
	grep -v "^>" $out.fna >> $out-concat.fna
	# compare contigs to ref
	nucmer -maxmatch -c 100 -p ${out}-concat $ref ${out}-concat.fna
	mummerplot --color -c -t png -p ${out}-concat-cov ${out}-concat.delta
	mummerplot --color    -t png -p ${out}-concat-dp  ${out}-concat.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-concat.delta > ${out}-concat.coords
	mapview -n 1 -f png -p ${out}-concat-map ${out}-concat.coords
	
	show-tiling ${out}-concat.delta > ${out}-concat.tiling
	mummerplot --color -c -t png -p ${out}-concat-tiling-cov ${out}-concat.tiling
	
	# compare ref to contigs
	nucmer -maxmatch -c 100 -p ${out}-tacnoc ${out}-concat.fna $ref
	mummerplot --color -c -t png -p ${out}-tacnoc-cov ${out}-tacnoc.delta
	mummerplot --color    -t png -p ${out}-tacnoc-dp  ${out}-tacnoc.delta

	show-coords -r -c -l -L 100 -I 50 ${out}-tacnoc.delta > ${out}-tacnoc.coords
	mapview -n 1 -f png -p ${out}-tacnoc ${out}-tacnoc.coords

	show-tiling ${out}-tacnoc.delta > ${out}-tacnoc.tiling
	mummerplot --color -c -t png -p ${out}-tacnoc-tiling-cov ${out}-tacnoc.tiling
    fi

    cd ..
    #set +x
}

function last_reference_align() {
    if [ "$LAST" != "yes" ] ; then return ; fi
    #set -x
    ref=`abspath $1`
    probe=`abspath $2`

    f="${probe##*/}"
    ext="${probe##*.}"
    out="${f%.*}"

    if [ ! -e $ref -o ! -e $probe ] ; then
        echo "Error, $ref and $probe must exist"
	return
    fi
    if [ ! -d last ] ; then
    	mkdir last
    fi
    echo "Entering last"
    cd last

    if [ ! -e reference.fna ] ; then
    	ln -sf `relpath $ref` reference.fna
    	#ln -sf $ref reference.fna
        #symlinks -c .
    fi
    if [ ! -e $out.fna ] ; then
        if [ ! -e $out.fastq ] ; then
    	    ln -sf `relpath $probe` $out.fastq
    	    #ln -sf $probe probe.fastq
	fi
        echo "creating $out.fna from $out.fastq"
        fq2f $out.fastq > $out.fna
    fi
    if [ ! -e refdb.prj ] ; then
    	echo "Making default db"
    	lastdb -c refdb reference.fna
    fi
    if [ ! -e $out.maf ] ; then
    	echo "aligning"
	# we use here the less stringent db and fasta without quality codes
	#	to get a minimum restriction (maximal information) alignment
    	lastal -Q0 refdb $out.fna > $out.maf
    	maf-convert html $out.maf > $out.html
    	maf-convert tab  $out.maf > $out.tab
    fi
    if [ ! -e $out.png ] ; then
        last-dotplot $out.maf $out.png
    fi
    
    # to require strong alignments:
    if [ ! -e refdbs.prj ] ; then
    	echo "Making strong db"
    	lastdb -m1111110 refdbs reference.fna
    fi
    if [ ! -e $out.strong.maf ] ; then
        # to indicate fastq-sanger format use -Q1, to indicate fasta -Q0
        echo "aligning with strong DB"
        lastal -Q1 refdbs $out.fastq > $out.strong.maf
    fi
    if [ ! -e $out.strong.png ] ; then
        last-dotplot $out.strong.maf $out.strong.png
    fi
    if [ ! -e $out.e1000.maf ] ; then
        # to ask for score > 1000
        echo "aligning score > 1000 with strong DB"
        lastal -Q1 -e1000 refdbs $out.fastq > $out.e1000.maf
    fi
    if [ ! -e $out.e1000.png ] ; then
        last-dotplot $out.e1000.maf $out.e1000.png
    fi
    if [ ! -e $out.1.maf ] ; then
	#set -x

        # to look for unique best alignments
        echo "selecting unique alignments"
	#
	# for some reason the fastq-based alignments fail to compute
        #last-split -g refdbs $out.strong.maf > $out.strong.1.maf
        last-split $out.strong.maf > $out.strong.1.maf
        last-dotplot $out.strong.1.maf $out.strong.1.png

        #last-split -g refdbs $out.e1000.maf > $out.e1000.1.maf
        last-split $out.e1000.maf > $out.e1000.1.maf
        last-dotplot $out.e1000.1.maf $out.e1000.1.png

	# but the fasta-based alignment has no problem 
	#last-split -g refdb  $out.maf > $out.1.maf
	last-split $out.maf > $out.1.maf
        last-dotplot $out.1.maf $out.1.png

        # sort by seqname, strand, start, end, of the top sequence
        maf-sort $out.1.maf > $out.1.srt.maf
	maf-convert tab $out.1.srt.maf > $out.1.srt.maf.tab
	
	#set +x
    fi
    for i in *.maf ; do
        echo -n $i
        if [ ! -e `basename $i maf`png ] ; then
            echo " ..plotting"# `pwd`/$i
            last-dotplot $i `basename $i maf`png
        else
            echo "."
        fi
    done
    
    cd ..
}

assemble_genome_from_paired_reads() {
    ref=$1
    r1=$2
    r2=$3

    if [ ! -e $ref -o ! -e $r1 -o ! -e $r2 ] ; then
        echo "Error: missing input file. Check that"
	echo "    $r"
	echo "    $r1"
	echo "    $r2"
	echo "do exist!"
	exit
    fi

    if [ "$BWA" = "yes" ] ; then
        echo "Doing BWA alignment"
        bwa_reference_assemble $ref $r1 $r2
	cd bwa
        echo "Analyzing SAM file"
        analyze_sam_file aligned.bwa.sam reference.fna
        echo "Aligning with MUMMER3"
        mummer3_align reference.fna aligned.bwa.cons.fastq
	echo "Aligning with NUCMER"
        nucmer_align reference.fna aligned.bwa.cons.fastq
	echo "Aligning with LAST"
        last_reference_align reference.fna aligned.bwa.cons.fastq
        cd ..
    fi
    if [ "$BOWTIE2" = "yes" ] ; then
        bowtie2_reference_assemble $ref $r1 $r2
	cd bowtie2
        analyze_sam_file aligned.bowtie.sam reference.fna
        mummer3_align reference.fna aligned.bowtie.cons.fastq
        nucmer_align reference.fna aligned.bowtie.cons.fastq
        last_reference_align reference.fna aligned.bowtie.cons.fastq
        cd ..
    fi
    if [ "$LASTASS" = "yes" ] ; then
        last_reference_assemble $ref $r1 $r2
	cd last
        analyze_sam_file paired.sam reference.fna
        mummer3_align reference.fna paired.cons.fastq
        nucmer_align reference.fna paired.cons.fastq
        last_reference_align reference.fna paired.cons.fastq
        cd ..
    fi
    if [ "$MAQ" = "yes" ] ; then
        maq_reference_assemble $ref $r1 $r2
        cd maq
        #mummer3_align reference.fna cns.fna
        #nucmer_align reference.fna cns.fq
        #last_reference_align reference.fna cns.fq
        cd ..
    fi
    if [ "$A5" = "yes" ] ; then
        a5_denovo_assemble $r1 $r2
        cd a5
	# WARNING!!! XXX JR XXX
	# this is TRICKY and depends STRONGLY on the current setup
        mummer3_align ../reference.fna genome.final.scaffolds.fastq
        nucmer_align ../reference.fna genome.final.scaffolds.fastq
        last_reference_align ../reference.fna genome.final.scaffolds.fastq
        cd ..
    fi
}
function usage()
{
    echo "Usage: $0 -r -i -b -B -M -A -m -l -n -v -h"
    echo ""
    echo "       $0 --roche454 --illumina --bwa --bowtie2 --maq --a5"
    echo "           --mummer3 --last --nucmer"
    echo ""
    echo "    -r --roche454    Reads are in Roche 454 format"
    echo "    -i --illumina    Reads are in illumina format"
    echo "                     These two are mutually exclusive"
    echo ""
    echo "    -b --bwa         Align mate pairs using BWA"
    echo "    -B --bowtie2     Align mates using Bowtie-2"
    echo "    -L --lastass     Align mates using the LAST aligner"
    echo "    -M --maq         Align mates using MAQ"
    echo "                     These three require a reference genome"
    echo "                     sequence in a subdirectory named 'refs'"
    echo "                     and in a file called 'reference.fna' in"
    echo "                     fasta format"
    echo "    -A --a5          Do a 'de novo' assemmbly using A5"
    echo "                     This option does not need a reference"
    echo "                     unless it is coupled with a comparison option"
    echo "    -m --mummer3     Compare obtained consensus against reference"
    echo "                     using MUMMER3"
    echo "    -l --last        Compare obtained consensus against reference"
    echo "                     using LAST aligner"
    echo "    -n --nucmer      Compare obtained consensus against reference"
    echo "                     using NUCMER"
    echo "                     The comparison options require a reference"
    echo "                     genome in fasta format, available as "
    echo "                     refs/reference.fna"
    echo "    -h --help        Pint this help and exit"
    echo "    -v --verbose     Output additional information"
    echo ""
    echo "    -S seqdir        A directory containing the sequenced reads"
    echo "                     in both directions, named *R[12]_001.fastq"
    echo "                     If none is specified, we'll use './orig'"
    echo "    -R ref.fna       The reference sequence to use for assembly"
    echo "                     and comparison"
    echo "                     If none is specified, we'll use './refs/reference.fna'"
    echo ""
    echo "Input reads must be in a subdirectory named '\$seqdir' and they"
    echo "    must conform to the naming convention *_R[12]_001.fastq"
    echo "    They can be optionally compressed with gzip"
    echo ""
    echo "Output will be saved in a directory named \$seqdir_x_\$ref"
    return
}

# --------------------------------------------------------------------------
#                              ASSEMBLE MATES
# --------------------------------------------------------------------------

# Parse the command line
# ----------------------
# Note that we use `"$@"' to let each command-line parameter expand to a 
# separate word. The quotes around `$@' are essential!
# We need TEMP as the `eval set --' would nuke the return value of getopt.
TEMP=`getopt -o hvribBLMAmlnS:1:2:R: \
     --long help,verbose,roche454,illumina,bwa,bowtie2,lastass,maq,a5,mummer3,last,nucmer,sequences:,R1:,R2:,reference: \
     -n "$0" -- "$@"`

# an invalid option was given
if [ $? != 0 ] ; then usage >&2 ; exit 1 ; fi

# Note the quotes around `$TEMP': they are essential!
eval set -- "$TEMP"

while true ; do
        case "$1" in
                -h|--help) 
                    usage ; return ; shift ;;

                -v|--verbose) 
                    VERBOSE=$(( $VERBOSE + 1 )) ; shift ;;
                    
                -r|--roche454)
                    SEQUENCE_IS_454="yes"
                    SEQUENCE_IS_ILLUMINA="no"
                    shift ;;
                
                -i|--illumina)
                    SEQUENCE_IS_454="no"
                    SEQUENCE_IS_ILLUMINA="yes"
                    shift ;;
                    
                -b|--bwa) 
                    BWA="yes" ; shift ;;
                
                -B|--bowtie2)
                    BOWTIE2="yes" ; shift ;;
                
                -L|--lastass)
                    LASTASS="yes" ; shift ;;
                
                -M|--maq)
                    MAQ="yes" ; shift ;;
                
                -A|--a5)
                    A5="yes" ; shift ;;
                
                -m|--mummer3)
                    MUMMER="yes" ; shift ;;
                
                -l|--last)
                    LAST="yes" ; shift ;;
                
                -n|--nucmer)
                    NUCMER="yes" ; shift ;;

		-S|--sequences)
                    SEQDIR=$2 ; shift 2 ;;
                
		-1|--R1)
                    R1=$2 ; shift 2 ;;
                
		-2|--R2)
                    R2=$2 ; shift 2 ;;
                
                -R|--reference)
                    REFFNA=$2 ;  shift 2 ;;
                    
                --) shift ; break ;;
                *) echo "Internal error!" >&2 ; usage ; exit 1 ;;
        esac
done

if [ $VERBOSE -gt 0 ] ; then
cat <<END
Selected values:

SEQUENCE_IS_454="$SEQUENCE_IS_454"		# -r --roche454
SEQUENCE_IS_ILLUMINA="$SEQUENCE_IS_ILLUMINA"	# -i --illumina
BWA="$BWA"			# -b --BWA
BOWTIE2="$BOWTIE2"		# -B --Bowtie2
MAQ="$MAQ"			# -M --MAQ
LASTASS="$LASTASS"		# -L --lastass
A5="$A5"			# -A --A5
MUMMER="$MUMMER"		# -m --mummer3
LAST="$LAST"			# -l --last
NUCMER="$NUCMER"		# -n --nucmer
VERBOSE="$VERBOSE"              # -v --verbose
SEQDIR="$SEQDIR"		# -S --sequences seqdir
REFFNA="$REFFNA"		# -R --reference ref.fna
R1="$R1"			# -1 --R1
R2="$R2"			# -2 --R2

END
fi

# Prepare input
# =============

# Set up reads files
# Absolutize R1 and R2 so we can move around
R1=`abspath $R1`
R2=`abspath $R2`
READS="orig"

if [ ! -d "$SEQDIR" ] ; then
    # if the dir does not exist
    if [ ! -e "$R1" -o ! -e "$R2" ] ; then
        echo "ERROR: We cannot find the reads in either"
	echo "    $SEQDIR"
	echo "    $R1" or "$R2"
        exit
    #else use $R1 and $R2
    fi
else
    READS="$SEQDIR"
    # if the dir exists, 
    if [ ! -e "$R1" -o ! -e "$R2" ] ; then
        # check if we can substitute for the missing $R1 and $R2
	if [ -e $SEQDIR/*_R1_001.fastq -a \
	     -e $SEQDIR/*_R2_001.fastq ] ; then
    	    R1="../$SEQDIR/*_R1_001.fastq"
	    R2="../$SEQDIR/*_R2_001.fastq"
	elif [ -e $SEQDIR/*_R1_001.fastq.gz -a \
	       -e $SEQDIR/*_R2_001.fastq.gz ] ; then
    	    R1="../$SEQDIR/*_R1_001.fastq.gz"
	    R2="../$SEQDIR/*_R2_001.fastq.gz"
	else
            echo "ERROR: We cannot find the reads in either"
	    echo "    $SEQDIR"
	    echo "    $R1" or "$R2"
            exit	
	fi
    #else use $R1 and $R2
    fi
fi

# find out the name of the reference (to be used later)
REFD=`dirname "$REFFNA"`
REFF="${REFFNA##*/}"
REF="${REFF%.*}"


# Check existence of reference sequence
REFERENCE=`abspath $REFFNA`
if [ ! -e $REFERENCE ] ; then
    echo "We need a reference sequence"
    echo "for reference assemblers and ulterior comparisons"
    echo "$REFFNA"
    echo "$REFERENCE"
    exit
fi


# Prepare FASTQ files
# ===================

if [ ! -d fastq ] ; then
    mkdir -p fastq
fi
cd fastq
echo "Populating fastq directory with initial reads"
if [ ! -e R1_001.fastq ] ; then
    if [ "${R1##*.}" == "gz" ] ; then
	zcat $R1 > R1_001.fastq
    else
	ln -s $R1 R1_001.fastq
    fi
fi
if [ ! -e R2_001.fastq ] ; then
    if [ "${R2##*.}" == "gz" ] ; then
	zcat $R2 > R2_001.fastq
    else
	ln -s $R2 R2_001.fastq
    fi
fi
cd ..


# Quality control
# ===============

if [ ! -d qc ] ; then
    mkdir -p qc
fi
cd qc
echo "Doing quality control"

# Check for mated pairs
#
if [ ! -e mates.log ] ; then
    ln -s ../fastq/R[12]_001.fastq .
    echo "Checking number of mates in each file (and in both)"
    cat R[12]_001.fastq | grep '^@.*:.*:.*:.*:.* .*:.*:.*:.*' \
        | cut -d' ' -f 1 | sort | uniq -d > mates.lst

    cat R1_001.fastq | grep '^@.*:.*:.*:.*:.* .*:.*:.*:.*' | wc -l > mates.log
    cat R2_001.fastq | grep '^@.*:.*:.*:.*:.* .*:.*:.*:.*' | wc -l >> mates.log
    wc -l mates.lst >> mates.log
    echo "" >> mates.log
    echo "The three counts should be the same" >> mates.log
    cat mates.log
fi

# FastX analysis of reads
#
if [ ! -e R1_001.qc ] ; then
    echo "Running FastX's fastqqc on R1"
    fastqqc R1_001.fastq > R1_001.qc
fi
if [ ! -e R2_001.qc ] ; then
    echo "Running FastX's fastqqc on R2"
    fastqqc R2_001.fastq > R2_001.qc
fi


# NGS_QC_Toolkit
#	Filter HQ reads
#
if [ "$SEQUENCE_IS_454" = "yes" ] ; then
    echo "Filtering using 454 QC"
    # R1
    if [ ! -d R1 ] ; then
        echo "NGS-QC filtering R1"
        perl ~/contrib/NGSQCToolkit_v2.3.3/Format-converter/FastqTo454.pl \
	        -i R1_001.fastq -o R1
        cd R1
        perl ~/contrib/NGSQCToolkit_v2.3.3/QC/454QC_PE.pl \
	        -i R1_001.fastq_fna  R1_001.fastq_qual 2
        cd ..
    fi

    # R2
    if [ ! -d R2 ] ; then
        echo "NGS-QC filtering R2"
        perl ~/contrib/NGSQCToolkit_v2.3.3/Format-converter/FastqTo454.pl \
	        -i R2_001.fastq -o R2
        cd R2
        perl ~/contrib/NGSQCToolkit_v2.3.3/QC/454QC_PE.pl \
	        -i R2_001.fastq_fna  R2_001.fastq_qual 2
        cd ..
    fi

    # all together now
    if [ ! -d R12 ] ; then
        echo "joining both read files"
        mkdir R12
        cat R1/R1_001.fastq_fna R2/R2_001.fastq_fna > R12/reads_R12_001.fastq_fna
        cat R1/R1_001.fastq_qual R2/R2_001.fastq_qual > R12/reads_R12_001.fastq_qual
        cd R12
        perl ~/contrib/NGSQCToolkit_v2.3.3/QC/454QC_PE.pl \
	        -i reads_R12_001.fastq_fna  reads_R12_001.fastq_qual 2
    fi

    # Use QIIME converter to generate fastq files
    #	we need to add -F -b to ensure trailing text after ' ' is included.
    if [ ! -d ../fastq/454_filtered ] ; then
        mkdir ../fastq/454_filtered
        echo "creating filtered FastQ files"
        echo "R1"
        if [ ! -e ../fastq/454_filtered/*R1_001.fastq ] ; then
            convert_fastaqual_to_fastq.py -f R1/454QC_Filtered_files/*R1_001.fastq_fna_filtered -q R1/454QC_Filtered_files/*R1_001.fastq_qual_filtered -F -b -o ../fastq/454_filtered
        fi
        echo "R2"
        if [ ! -e ../fastq/454_filtered/*R2_001.fastq ] ; then
            convert_fastaqual_to_fastq.py -f R2/454QC_Filtered_files/*R2_001.fastq_fna_filtered -q R2/454QC_Filtered_files/*R2_001.fastq_qual_filtered -F -b -o ../fastq/454_filtered
        fi
        echo "R12"
        if [ ! -e ../fastq/454_filtered/*_R12_001.fastq ] ; then
            convert_fastaqual_to_fastq.py -f R12/454QC_Filtered_files/*R12_001.fastq_fna_filtered -q R12/454QC_Filtered_files/*R12_001.fastq_qual_filtered -F -b -o ../fastq/454_filtered
        fi
    fi
fi	# SEQUENCE_IS_454 = yes
if [ "$SEQUENCE_IS_ILLUMINA" = "yes" ] ; then
    echo "Filtering using Illumina QC"
    # Filter using Illumina filter
    if [ ! -d  IlluQC_Filtered_files ] ; then
        perl ~/contrib/NGSQCToolkit_v2.3.3/QC/IlluQC_PRLL.pl \
	    -pe R1_001.fastq R2_001.fastq 2 a -c 8
    else
    	echo "An Illumina quality directory already exists"
    fi
    if [ ! -d ../fastq/Ill_filtered ] ; then
        mkdir ../fastq/Ill_filtered 
    fi
    cd ../fastq/Ill_filtered
    for i in ../../qc/IlluQC_Filtered_files/*.fastq_filtered ; do
        if [ ! -h `basename $i _filtered` ] ; then
            ln -s $i `basename $i _filtered`
        fi
    done
    cd ../../qc
fi

# exit qc directory
cd ..

# Analyze sequences
#
#	with R1.fastq, R2.fastq and reference.fna
#
if [ "$SEQUENCE_IS_454" = "yes" ] ; then
    echo "Analyzing 454 sequences"
    if  [ ! -d "${READS}_x_${REF}_454" ] ; then
        mkdir "${READS}_x_${REF}_454"
    fi
    REFERENCE=`abspath $REFFNA`
    cd "${READS}_x_${REF}_454"
    if [ ! -e reference.fna ] ; then
        if [ -e $REFERENCE ] ; then
	    ln -s `relpath $REFERENCE` reference.fna
        else
            echo "We need a reference sequence"
            echo "for reference assemblers and ulterior comparisons"
	    exit
        fi
    fi
    if [ ! -e R1.fastq ] ; then
        if [ -e ../fastq/454_filtered/*R1_001.fastq ] ; then
            ln -s ../fastq/454_filtered/*R1_001.fastq R1.fastq
        else
            echo "We need a filtered reads file to run 454 assemblers"
            exit
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        if [ -e ../fastq/454_filtered/*R2_001.fastq ] ; then
            ln -s ../fastq/454_filtered/*R2_001.fastq R2.fastq
        else
            echo "We need a filtered reads file to run 454 assemblers"
            exit
        fi
    fi
fi	# SEQUENCE_IS_454 = yes

if [ "$SEQUENCE_IS_ILLUMINA" = "yes" ] ; then
    echo "Analyzing Illumina sequences"
    if  [ ! -d "${READS}_x_${REF}" ] ; then
        mkdir "${READS}_x_${REF}"
    fi
    REFERENCE=`abspath $REFFNA`
    cd "${READS}_x_${REF}"
    if [ ! -e reference.fna ] ; then
        if [ -e $REFERENCE ] ; then
	    ln -s `relpath $REFERENCE` reference.fna
        else
            echo "We need a reference sequence"
            echo "for reference assemblers and ulterior comparisons"
	    exit
        fi
    fi
    echo "Preparing FastQ input files..."
    if [ ! -e R1.fastq ] ; then
        #if [ -e ../fastq/Ill_filtered/*_R1_001.fastq ] ; then
        if [ "$(ls -A ../fastq/Ill_filtered/R1_001.fastq)" ] ; then
            ln -s ../fastq/Ill_filtered/R1_001.fastq R1.fastq
            # checking if this works for multiple files
            #cat ../fastq/Ill_filtered/R1_001.fastq > R1.fastq
        else
            echo "We need a filtered reads file (1)"
            exit
        fi
    fi
    if [ ! -e R2.fastq ] ; then
        #if [ -e ../fastq/Ill_filtered/*R2_001.fastq ] ; then
        if [ "$(ls -A ../fastq/Ill_filtered/R2_001.fastq)" ] ; then
            ln -s ../fastq/Ill_filtered/R2_001.fastq R2.fastq
            #cat ../fastq/Ill_filtered/R2_001.fastq > R2.fastq
        else
            echo "We need a filtered reads file (2)"
            exit
        fi
    fi
fi


assemble_genome_from_paired_reads reference.fna R1.fastq R2.fastq

cd ..
exit

