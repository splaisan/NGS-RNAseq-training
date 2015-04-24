#!/bin/bash
# cutadapt_both.sh
# apply to SRR1039509_1/2

# takes a fastQ file from $1 and performs 'cutadapt' cleaning
if [ $# -ge 1 ]
then
echo "Usage: ${0##*/} <fastQ file_1> (second fastQ is deduced)"
exit
fi

## run using gnu-parallel with 5 concurrent jobs like
## find SRP012167_fastq/*_1.fastq.gz \
##	| parallel --no-notice \
##		--workdir . --tmpdir ./tmp \
##		-j 5 scripts/cutadapt_both.sh {}

outfolder="trimmed-SRR1039509_fastq"
mkdir -p ${outfolder}

# save tmp files elsewhere
tmpdir=./tmp
mkdir -p ${tmpdir}

# process the FastQ pair
fq1=${1:-"SRR1039509_fastq/SRR1039509_1.fastq.gz"}
name1=$(basename ${fq1} .fastq.gz)
fq2=${fq1%%_1.fastq.gz}_2.fastq.gz
name2=$(basename ${fq2} .fastq.gz)

#echo "# decompressing input files"
#gzip -dc -f ${fq1} > ${tmpdir}/${name1}.fastq
#gzip -dc -f ${fq2} > ${tmpdir}/${name2}.fastq

echo "# removing adaptors from read pair: ${name1} and ${name2}"

# TruSeq Adapter, Index 1 (100% over 50bp)
# ADAPTER_FWD="GATCGGAAGAGCACACGTCTGAACTCCAGTCACATCACGATCTCGTATGC"

# TruSeq Adapter, Index 4 (100% over 51bp) - SRR1039509_1.fastq.gz
ADAPTER_FWD="ACACGTCTGAACTCCAGTCACTGACCAATCTCGTATGCCGTCTTCTGCTTGAAAAAAAAAAAA"
cutadapt -q 10 \
	-a ${ADAPTER_FWD} \
	--minimum-length 20 \
	--paired-output ${tmpdir}/${name2}.tmp.fastq \
	-o ${tmpdir}/${name1}.tmp.fastq \
	<(gzip -dc -f ${fq1}) \
	<(gzip -dc -f ${fq2}) \
	2>&1 | tee -a ${outfolder}/cutadapt-${name1%%_1}_stats.txt

echo "# done for ${name1} and forward adaptor"

# Illumina Single End PCR Primer 1 (100% over 50bp)
# ADAPTER_REV="GATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCG"

# Illumina Single End PCR Primer 1 (100% over 45bp) - SRR1039509_2.fastq.gz
ADAPTER_REV="GTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAAAAAAAAAAAAAA"
cutadapt -q 15 \
	-a ${ADAPTER_REV} \
	--minimum-length 20 \
	--paired-output ${tmpdir}/${name1}.trimmed.fastq \
	-o ${tmpdir}/${name2}.trimmed.fastq \
	${tmpdir}/${name2}.tmp.fastq \
	${tmpdir}/${name1}.tmp.fastq \
	2>&1 | tee -a ${outfolder}/cutadapt-${name1%%_1}_stats.txt

echo "# done for ${name2} and reverse adaptor"

echo "# compressing results and removing tmp files"

bgzip -c ${tmpdir}/${name1}.trimmed.fastq > \
	${outfolder}/${name1}.fastq.gz \
	&& rm ${tmpdir}/${name1}.trimmed.fastq

bgzip -c ${tmpdir}/${name2}.trimmed.fastq > \
	${outfolder}/${name2}.fastq.gz \
	&& rm ${tmpdir}/${name2}.trimmed.fastq

echo "# done all"
