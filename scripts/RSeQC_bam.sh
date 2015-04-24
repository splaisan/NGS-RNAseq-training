#!/bin/bash
# RSeQC_bam.sh

# takes a BAM file from $1 and performs several RSeQC analyses
if [ $# -lt 1 ]; then
	echo "Usage: ${0##*/} <BAM file>"
	exit
fi

# process only one BAM file
# batch can be done under bash with a for loop
bam=$1

# help: batch apply this script to all bams with
## for b in markdup_bams/*.bam; do \
##   echo "# analyzing ${b}" \
##   scripts/RSeQC_bam.sh $b\
##   done

qc_results="RSeQC_results"
mkdir -p ${qc_results}

echo "# analysing: ${bam}"
name=$(basename "${bam%%-*.bam}")

### DOWNLOAD and expand reference transcript models ###
# get reference files from the RSeQC server if absent locally
mkdir -p ref

# refseq models
if [[ ! -f ref/hg19_RefSeq.bed.gz ]]; then
	wget -P ref \
	http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/hg19_RefSeq.bed.gz
fi

# ensembl models
if [[ ! -f ref/hg19_Ensembl.bed.gz ]]; then
	wget -P ref \
	http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/hg19_Ensembl.bed.gz
fi

# we choose the EnsEMBL reference model
orirefgene="hg19_Ensembl.bed.gz"
refgene="hg19_ensGene.bed"

# expand ensGene and sort chromosome records in ascending order
# and create chr22 subset
gzip -cd ref/${orirefgene} | \
	sort -k 1V,1 -k2n,2 -k 3n,3 > ref/${refgene} &&
	grep "^chr22" ref/${refgene} > ref/chr22-${refgene}

##############################
# perform selected RSeQC tests
# more tools exist, please refer to the online documentation.

# select which reference transcriptome to use
# or chr22 subset
selref=ref/${refgene}
#selref=ref/chr22-${refgene}

# prefix to save results
res=${qc_results}/${name}

# build the 'infer_experiment' command
# speculate how RNA-seq sequencing were configured
cmd="infer_experiment.py \
	-r ${selref} \
	-i ${bam} \
	-s 200000 \
	> ${res}_inferred.txt \
	2>&1 | tee -a ${qc_results}/RSeQC-${name}_log.txt"

echo "# ${cmd}"
eval ${cmd}

# build the 'inner_distance.py' command
# calculate the inner distance (or insert size) between two paired RNA reads
# The inner_distance might be a negative value if paired reads do overlap.
cmd="inner_distance.py \
	-r ${selref} \
	-i ${bam} \
	-o ${res}_inner_distance \
	-k 1000000 \
	-l -250 \
	-u 250 \
	-s 5 \
	-q 30 \
	2>&1 | tee -a ${qc_results}/RSeQC-${name}_log.txt"

echo "# ${cmd}"
eval ${cmd}

# build the 'junction_saturation' command
# check if sequencing depth is deep enough for alternative splicing analyses.
cmd="junction_saturation.py \
	-r ${selref} \
	-i ${bam} \
	-o ${res} \
	-l 5 \
	-u 100 \
	-m 50 \
	-v 1"

echo "# ${cmd}"
eval ${cmd}

# build the 'geneBody_coverage' command
# compute coverage for all transcript models in percentiles
cmd="geneBody_coverage.py \
	-r ${selref} \
	-i ${bam} \
	-o ${res} \
	2>&1 | tee -a ${qc_results}/RSeQC-${name}_log.txt"

echo "# ${cmd}"
eval ${cmd}

# determine reads duplication rate with 'read_duplication.py'
cmd="read_duplication.py \
	-i ${bam} \
	-o ${res} \
	-u 500 \
	2>&1 | tee -a ${qc_results}/RSeQC-${name}_log.txt"

echo "# ${cmd}"
eval ${cmd}

echo "# end of RseQC analysis"

# rem: the following 'ImageMagick' command converts all PDF outputs to PNG
# and replaces internal '.' that would disturb some packages by '_'
# not run:
# for i in *.pdf; do name=${i/./_}; convert "$i" "${name%.pdf}".png; done

