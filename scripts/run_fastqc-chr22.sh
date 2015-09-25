#!/bin/bash
# run_fastqc.sh

# chr22 data
# we read and save data in the $BASE folder
echo $BASE

# first line processes all fastq while second takes only one dataset
#export RAWDATA=$BASE/SRP033351-chr22_fastq
export RAWDATA=$BASE/SRR1039509-chr22_fastq

export READQC=$BASE/SRR1039509-chr22_read_qc-Results

mkdir -p $READQC

# create an empty error log
cat /dev/null > $READQC/error.log

# max threads
thr=2

# process all fastq in folder
for fq in $RAWDATA/*.fastq.gz; do
	name=${fq##*/}
	prefix=${name%%.fastq.gz}

	# perform a full QC control on the reads (-q for quiet, 4 threads)
	# then convert results to PDF
	echo "# run fastqc on ${fq}"
	cmd="fastqc -f fastq --noextract -t ${thr} -o $READQC -q ${fq}"

	echo "# ${cmd}"
	eval ${cmd}

	## add conversion to PDF
	# htmldoc --webpage -f ${prefix}_fastqc.pdf ${prefix}_fastqc.html

	## Use now another package 'FASTX toolkit' to produce additional QC data
	# generate statistics for plots
	# secret option '-Q 33' was added to fit with the phred scale
	echo "# checking: ${name}"
	zcat ${fq} | fastx_quality_stats -Q 33 -o $READQC/${name}_stats.txt \
		2>> $READQC/error.log

	# plot from the text summary file
	fastq_quality_boxplot_graph.sh \
		-i $READQC/${name}_stats.txt \
		-o $READQC/${prefix}_boxplot.png

	fastx_nucleotide_distribution_graph.sh \
		-i $READQC/${name}_stats.txt \
		-o $READQC/${prefix}_nuclplot.png

	# also plot normalized base frequency using R (script code in appendix)
	scripts/avgQdist2linePlot.R $READQC/${name}_stats.txt $READQC

done

# convert results to a nice PDF file
# requires htmldoc and dependencies installed on your machine
# htmldoc --webpage \
#     --browserwidth 800 \
#     --fontsize 7 \
#     -f ${outfolder}/fastqc_report.pdf \
#     ${outfolder}/fastqc_report.html10
