
targets="align/*.recalibrated.bam"

data=$HOME/work/lorena/data/hg19
scripts=$HOME/work/lorena/scripts
out=cnv.cnvkit

baits_bed=chrexome.intervals.bed

if [ ! -d $out ] ; then
    mkdir -p $out
fi

ln -s $data/refFlat.txt .
ln -s $data/ucsc.hg19.fasta.* .
ln -s $scripts/cnvkit.py .

./cnvkit.py batch \
	align/*.recalibrated.bam \
	-n \
	-t $baits_bed \
	-f $data/ucsc.hg19.fasta \
	--annotate refFlat.txt \
	--output-reference $out/reference.cnn \
	--output-dir $out \
	-p 0 \
	--scatter \
	--diagram \
	-d $out

