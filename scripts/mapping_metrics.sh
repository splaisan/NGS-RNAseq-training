#! /usr/bin/env bash
## script: 'mapping_metrics.sh'
## ©SP-BITS, 2015-01-05

# required:
# samtools 1.2
# Picard Version (broad build): 1.129

maxmem=4G

infile=$1
filename=$(basename ${infile})
outbase="${infile:0:${#infile} - ${#filename}}"

## samtools flagstats ##
cmd="samtools flagstat \
	${infile} \
	>${infile%%.bam}_flagstat.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools stats ##
cmd="samtools stats \
	${infile} \
	>${infile%%.bam}_stats.txt"

echo
echo "# ${cmd}"
eval ${cmd}

## samtools plot-bamstats ##
cmd="plot-bamstats \
	${infile%%.bam}_stats.txt \
	-p ${infile%%.bam}_plots/"

echo
echo "# ${cmd}"
eval ${cmd}

## plot mapping quality score distribution
cmd="java -Xmx${maxmem} -jar $PICARD/picard.jar \
	QualityScoreDistribution \
	CHART=${infile%%.bam}-QS_distribution.pdf \
	I=${infile} \
	O=${infile%%.bam}-mappings_QS.txt"

echo
echo "# ${cmd}"
eval ${cmd}
