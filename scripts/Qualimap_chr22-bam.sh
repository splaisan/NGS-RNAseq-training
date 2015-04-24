#!/bin/bash
# Qualimap_bam.sh
# bam needs to be sorted by coordinate

# takes a BAM file from $1 and performs several qualimap analyses
if [ $# -lt 1 ]; then
	echo "Usage: ${0##*/} <BAM file>"
	exit
fi

# process only one BAM file
# batch can be done under bash with a for loop
bam=$1
echo "# analysing: ${bam}"
name=$(basename "${bam%.bam}")

# create folders
qc_reports="Qualimap_results/chr22-reports/${name}"
mkdir -p ${qc_reports}

qc_counts="Qualimap_results/chr22-counts"
mkdir -p ${qc_counts}

# we choose the EnsEMBL reference model
#gtffile=$BOWTIE2_INDEXES/hg19_ensGene.gtf

# but we also have a chr22 subset file
gtffile=$BOWTIE2_INDEXES/chr22-hg19_ensGene.gtf

##############################
# perform selected Qualimap tests
# more tools exist, please refer to the online documentation.

# default names
sbam=${bam}
qbam=${bam}

# test if bam is query or coordinate-sorted
order="$(samtools view -H ${bam} | head -1 | awk '{split($3, ord, ":"); print ord[2]}')"

if [ "$order" = "coordinate" ]; then

	echo "# sorting in queryname order"
	# redefine qbam
	qbam=${bam//.bam/-q.bam}
	cmd="java -jar $PICARD/picard.jar SortSam \
		I=${bam} \
		O=${qbam} \
		SO=queryname"
	eval ${cmd}

else

	echo "# sorting in coordinate order"
	# redefine sbam
	sbam=${bam//.bam/-s.bam}
	cmd="java -jar $PICARD/picard.jar SortSam \
		I=${bam} \
		O=${sbam} \
		SO=coordinate"
	eval ${cmd}

fi

#echo "# ${cmd}"
#eval ${cmd}

# run the BAMQC tool on coordinate-sorted data
cmd="$QUALIMAP/qualimap bamqc \
 	-p non-strand-specific \
 	-bam ${sbam} \
 	-c \
 	-gd HUMAN \
 	-gff ${gtffile} \
 	-hm 3 \
 	-nr 1000 \
 	-nt 4 \
 	-nw 400 \
 	-outdir ${qc_reports} \
 	-outformat PDF \
 	-outfile ${name}-Qualimap-bamqc-report.pdf \
	2>&1 | tee -a ${qc_reports}/qualimap-bamqc_log.txt"

echo "# ${cmd}"
eval ${cmd}

# print separator
A=$(printf "%50s\n")
echo ${A// /#}
echo

# run the RNASEQ tool on queryname-sorted data
cmd="$QUALIMAP/qualimap rnaseq \
 	-p non-strand-specific \
 	-pe \
 	-bam ${qbam} \
 	-a uniquely-mapped-reads \
 	-gtf ${gtffile} \
  	-outdir ${qc_reports} \
 	-oc ${qc_counts}/${name}_counts.txt \
 	-outformat PDF \
 	-outfile ${name}-Qualimap-rnaseq-report.pdf \
	2>&1 | tee -a ${qc_reports}/qualimap-rnaseq_log.txt"

echo "# ${cmd}"
eval ${cmd}

# deleted leftover sorted BAM
if [ -f ${bam//.bam/-s.bam} ]; then
	rm ${bam//.bam/-s.bam}
fi

if [ -f ${bam//.bam/-s.bam} ]; then
	rm ${bam//.bam/-q.bam}
fi

echo "# end of Qualimap analysis"
