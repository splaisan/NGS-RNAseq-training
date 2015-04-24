#!/bin/bash
# tophat_map-all.sh

# the script should be called from the base folder
# the base folder includes the 'reads' folder named below

## one sample chr22 subset
reads="SRR1039509-chr22_fastq"

results="tophat2.0.13_results_gtf"
mkdir -p ${results}

# folder and file options
fasta=$BOWTIE2_INDEXES/hg19.fa
index=$BOWTIE2_INDEXES/hg19
gtffile=$BOWTIE2_INDEXES/hg19_ensGene.gtf
gtfindex=$BOWTIE2_INDEXES/hg19_ensGene

# how hard can your computer work? (max cpu# that may be used)
nthr=1 # using 'nthr' processors
maxram="4G"

# on a stronger, this will run faster
# nthr=24 # using 'nthr' processors
# maxram="64G"

# map single ends in (un-guided|guided) modes (not|using) the gtf file
for fq1 in $(ls ${reads}/*_1.fastq.gz); do

 # get the base-name without end_suffix
 name=$(basename "${fq1}" "_1.fastq.gz")
 echo
 echo "# mapping the ${name} pair against the whole genome"

 # deduce second read file name
 fq2=$(echo "$fq1" | sed -r 's/_1.fastq.gz/_2.fastq.gz/')

 # guided mode with 'ensGene' as model
 # execute the command and write summary results in a 'log' file
 # first time -G ${gtffile} then --transcriptome-index ${gtfindex}
 cmd="tophat \
 	-p ${nthr} \
 	-o ${results}/${name}_gtf-mappings \
 	--transcriptome-index ${gtfindex} \
 	--no-coverage-search \
 	${index} \
 	${fq1} ${fq2} \
 	2>&1 | tee -a ${results}/tophat-gtf-log_${name}.txt"

 echo "# ${cmd}"
 eval ${cmd}

 # if tophat was a success, sort, rename, and index mapping results
 # run in background with good resource while starting the next loop on 1 thread
 if [ $? = 0 ]; then
  # use 4 threads and up to 4GB RAM (default are 1 thread and 768M)
  samtools sort -@ ${nthr} -m ${maxram} \
  	${results}/${name}_gtf-mappings/accepted_hits.bam \
  	${results}/${name}_gtf-tophat \
	&& samtools index ${results}/${name}_gtf-tophat.bam

  # also perform basic BAM QC on mappings
  # the code is provided in a separate script 'mapping_metrics.sh'
  scripts/mapping_metrics.sh ${results}/${name}_gtf-tophat.bam
 fi

 echo

 # print separator
 A=$(printf "%50s\n")
 echo ${A// /#}
 echo

done
