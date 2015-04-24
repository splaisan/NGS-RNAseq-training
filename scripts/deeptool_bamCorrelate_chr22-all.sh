#!/bin/bash
# deeptool_bamCorrelate_chr22-all.sh

# global correlation
# we now compare all BAM files at once to look for clustering between replicates
# we extract only chr22 mappings and compare to ensGene loci

# required:
# deeptool bamCorrelate all 16 files (chr22 only and from 'all' mappings)
# https://github.com/fidelram/deepTools/wiki/QC#bamCorrelate
# https://github.com/fidelram/deepTools/wiki/All-command-line-options

# speedup with --numberOfProcessors nthr matching your cpu#
# nthr=24
nthr=2

bamfolder=SRP033351-chr22_bam
results=bamCorrelate_chr22-all
mkdir -p ${results}

### DOWNLOAD and expand reference transcript models ###
# get reference files from the RSeQC server if absent locally
mkdir -p ref
# refseq models
if [[ ! -f ref/hg19_RefSeq.bed.gz ]]; then
 wget -P ref \
	http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/hg19_RefSeq.bed.gz \
	&& zgrep "^chr22" ref/hg19_RefSeq.bed.gz > ref/chr22-hg19_RefSeq.bed
fi
# ensembl models
if [[ ! -f ref/hg19_Ensembl.bed.gz ]]; then
 wget -P ref \
	http://dldcc-web.brc.bcm.edu/lilab/liguow/RSeQC/dat/hg19_Ensembl.bed.gz \
	&& zgrep "^chr22" ref/hg19_Ensembl.bed.gz > ref/chr22-hg19_Ensembl.bed
fi

# we choose the Ensembl standard
bedfile="ref/chr22-hg19_Ensembl.bed"

# create bam file list
bamfiles=$(ls ${bamfolder}/*.bam)
# deduce labels from files
l1=${bamfiles//${bamfolder}\/chr22_/}
labels=${l1//_all.bam/}

# build the deeptool command
corM=spearman # spearman or pearson

cmd="bamCorrelate bins \
	--bamfiles ${bamfiles} \
	--labels ${labels} \
	--fragmentLength 200 \
	--numberOfProcessors ${nthr} \
	--outFileCorMatrix ${results}/correlation_matrix.txt \
	--binSize 1000 \
	--corMethod ${corM} \
	-o ${results}/correlation_${corM}.pdf \
	2>&1 | tee -a ${results}/bamCorrelate_log.txt"

echo "# ${cmd}"
eval ${cmd}
