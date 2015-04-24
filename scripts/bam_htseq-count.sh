#!/bin/bash
# bam_htseq-count.sh

# run HTSeq on single BAM sample
# count only uniquely mapped mappings & alignQ > 10
# requires HTSeq (0.6.0+)
## requires sam data sorted by 'queryname' as in our markdup set

# help: batch apply this script to all bams with
## for b in tophat2.0.13_results/markdup_bams/*_all-tophat-q.bam; do \
##   echo "# analyzing ${b}" \
##   scripts/bam_htseq-count.sh $b\
##   done

# or even better run 5 jobs in 'parallel' if you have power
# find tophat2.0.13_results/markdup_bams/ -name "*_all_mdup.bam" | \
#	parallel --no-notice --workdir . --tmpdir ./tmp -j 5 \
#	scripts/bam_htseq-count.sh {}

# takes a single bam file from $1
if [ $# -lt 1 ]
then
echo "Usage: ${0##*/} <bam file>
	eg: tophat2.0.13-SRR1039509_results/markdup_bams/SRR1039509_all_mdup.bam"
exit
fi
# argument#1 gives the bam file
bam=$1
name=$(basename ${bam} ".bam")

# test if bam is queryname-sorted
order="$(samtools view -H ${bam} | head -1 | awk '{split($3, ord, ":"); print ord[2]}')"
echo "# order is $order"
echo

if [ "$order" = "coordinate" ]; then
	echo "# sorting in queryname order"
	# redefine qbam
	qbam=${bam//.bam/-q.bam}
	cmd="java -jar $PICARD/picard.jar SortSam \
		I=${bam} \
		O=${bam//.bam/-q.bam} \
		SO=queryname"
	eval ${cmd} && bam=${bam//.bam/-q.bam}
fi

# create output folder
outfolder="htseq_counts"
mkdir -p ${outfolder}

# select transcript annotation model
gtffile=$BOWTIE2_INDEXES/hg19_ensGene.gtf
# other options would be
#gtffile=$BOWTIE2_INDEXES/chr22-hg19_ensGene.gtf
#gtffile=$BOWTIE2_INDEXES/hg19_refGene.gtf

# decompress bam with 'samtools' and use default parameters
# possible IDs are:
## transcript_name (gene symbol)
## transcript_id (ENST-ID)
## exon_id (ENSE-ID)
## gene_id (ENSG-ID:default) ++

cmd="samtools view ${bam} | \
  	htseq-count \
  	-m intersection-strict \
  	-s no \
  	-a 10 \
  	-t exon \
  	-i gene_id \
  	- \
  	${gtffile} \
  	> ${outfolder}/${name}_counts.txt \
 	2> >(tee -a ${outfolder}/${name}_errorlog.txt >&2)"

echo "# ${cmd}"
eval ${cmd}
echo

# compute simple stats on ENSG transcripts with coverage
# use non-null data from col#2 when col#1 contains a valid "ENSG ID"
echo "# simple coverage stats for transcripts with read counts"
echo "# ${bam}"
awk 'BEGIN{FS="\t"; OFS="\t"}{if($1~/ENSG/ && $2>0) print $2}' \
	${outfolder}/${name}_counts.txt | qstats -s

