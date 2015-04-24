#!/bin/bash
# bam_markdups.sh

# mark optical and PCR duplicates from 'SO:coordinate'-sorted BAM file(s)
# re-sort output by read name ('SO:queryname')
# requires picard

# allocate as much RAM as your computer can
maxram=4G

# mark optical and PCR duplicates from position-sorted BAM file(s)
infolder=${1:-"tophat2.0.13_results"}
outfolder=${infolder}/markdup_bams
mkdir -p ${outfolder}

for b in ${infolder}/*_???-tophat.bam; do

	echo "# processing ${b}"
	# get prefix
	name=$(basename ${b} "-tophat.bam")

	# mark duplicates and index result
	cmd="java -Xmx${maxram} -jar $PICARD/picard.jar MarkDuplicates \
		I=${b} \
		O=${outfolder}/${name}_mdup.bam \
		M=${outfolder}/${name}_mdup.txt \
		AS=TRUE \
		2>&1 | tee -a ${outfolder}/MarkDuplicates_${name}_log.txt"

	echo "# ${cmd}"
	eval ${cmd}

	# re-sort by querynames to comply with HTSeq-count
	# if last command did exit without error ($? -eq 0)
	if [ $? -eq 0 ]; then
		cmd="java -Xmx${maxram} -jar $PICARD/picard.jar SortSam \
			I=${outfolder}/${name}_mdup.bam \
			O=${outfolder}/s${name}_mdup.bam \
			SO=queryname \
			CREATE_INDEX=TRUE \
			VALIDATION_STRINGENCY=LENIENT \
			2>&1 | tee -a ${outfolder}/MarkDuplicates_${name}_log.txt"

		echo "# ${cmd}"
		# delete original file if picard exists successfully
		eval ${cmd} && mv ${outfolder}/s${name}_mdup.bam \
		${outfolder}/${name}_mdup.bam

	else
		echo "# something went wrong!"
		exit 1
	fi

	# output results to screen from text file
	echo "# MarkDuplicates results for ${name}"
	head ${outfolder}/${name}_mdup.txt | \
		grep -v "^#" | transpose -t | column -t | \
		tee -a > ${outfolder}/${name}_mdup-counts.txt
	
	echo

done
