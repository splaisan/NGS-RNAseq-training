#!/bin/bash
# extract_chr22_to_fastq.sh

# extract chr22 bam records to new bam
# extract chr22 reads to new fastq files

# takes a BAM file from $1 and extracts paired chr22 reads to fastq's
if [ $# -lt 1 ]
then
echo "Usage: ${0##*/} <BAM file>"
exit
fi

# edit here
maxram="2G"

# get BAM from argument#1
bam=$1

# set names
name=$(basename ${bam} ".bam")

outbam=chr22_${name}
mkdir -p ${outbam}

outfastq=chr22_${name}_fastq
mkdir -p ${outfastq}

echo "# processing ${bam}"

# add suffix
sbam=${bam//.bam/_s.bam}

# sorted version required for region extraction
cmd="java -Xmx${maxram} -jar $PICARD/picard.jar SortSam \
	I=${bam} \
	O=${sbam} \
	SO=coordinate"

echo "# ${cmd}"
eval ${cmd}

# index the obtained data before extraction
# same as: samtools index ${sbam}
cmd="java -Xmx${maxram} -jar $PICARD/picard.jar BuildBamIndex \
	I=${sbam} \
	O=${sbam}.bai"

echo "# ${cmd}"
eval ${cmd}

# extract chr22 mappings and index results
# only mapped and properly paired reads (tag=0x2)
cmd="samtools view -b -f 0x2 ${sbam} chr22 \
	> ${outbam}/${outbam}.bam && \
	samtools index ${outbam}/${outbam}.bam"

echo "# ${cmd}"
eval ${cmd}

# re-sort by read names to group reads in pairs
cmd="java -Xmx${maxram} -jar $PICARD/picard.jar SortSam \
	I=${outbam}/${outbam}.bam \
	O=${outbam}/${outbam}-n.bam \
	SO=queryname"

echo "# ${cmd}"
eval ${cmd}

# convert back to fastq using picard
# be LENIENT for unpaired reads
cmd="java -Xmx${maxram} -jar $PICARD/picard.jar SamToFastq \
	I=${outbam}/${outbam}-n.bam \
	F=>(bgzip -c > ${outfastq}/${outfastq}_1.fastq.gz) \
	F2=>(bgzip -c > ${outfastq}/${outfastq}_2.fastq.gz) \
	RE_REVERSE=TRUE \
	VALIDATION_STRINGENCY=LENIENT"

echo "# ${cmd}"
eval ${cmd}

# delete sorted files to save space
rm ${sbam} ${sbam}.bai
rm ${outbam}/${outbam}-n.bam
			
