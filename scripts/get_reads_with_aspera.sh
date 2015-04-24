#!/bin/bash
# get_reads_with_aspera.sh

# required: SRAtoolkit installed and in PATH
# required: a functional aspera connect installation

## Aspera main parameters
# ÐQ (for adaptive flow control) Ð needed for disk throttling!
# ÐT to disable encryption
# Ðk1 enable resume of failed transfers
# Ðl (maximum bandwidth of request, try 200M and go up from there)
# Ðr recursive copy
# Ði <private key file>

## set the following two variables according to your own Aspera location
# this should match sra-connect installed for a single user under unix
# /home/<user>/.aspera/connect
exefile="~/.aspera/connect/bin/ascp"
sshcert="~/.aspera/connect/etc/asperaweb_id_dsa.openssh"
baseurl="anonftp@ftp-private.ncbi.nlm.nih.gov:"

PRJ=SRP033351
BIOPJ=PRJNA229998
# create a LIST of 16 files to download and proceed with aspera
LIST="SRR1039508 SRR1039509 SRR1039510 SRR1039511 SRR1039512 SRR1039513 \
SRR1039514 SRR1039515 SRR1039516 SRR1039517 SRR1039518 SRR1039519 \
SRR1039520 SRR1039521 SRR1039522 SRR1039523"

# adapt the following line if you need other reads
# it should end with the first letters of your SRA files

# create container for data and move into it
mkdir -p ${PRJ}_fastq/sra_downloads && cd ${PRJ}_fastq

# loop in the LIST and download one file at a time
for sra in $LIST; do

	# convert ID to URL using SRAtoolkit 'srapath'
	uri=$(srapath ${sra} | sed 's|http://ftp-trace.ncbi.nlm.nih.gov|baseurl|')

	## download SRA data
	cmd="${exefile} \
		-i ${sshcert} \
		-k 1 \
		-QTr \
		-l 10000m \
		${uri} \
		sra_downloads"

	echo "# $cmd"
	eval $cmd

	# test if transfer succeeded or die
	RESULT=$?

	if [ $RESULT -eq 0 ]; then
		echo "# Aspera transfer succeeded for ${sra}"
	else
		echo "# Aspera transfer failed for ${sra}, aborting!"
		exit 1
	fi

	## convert SRA to fastq
	# archive data in gzip format to save space
	fastq-dump --split-3 --gzip sra_downloads/${sra}.sra

	# test if conversion succeeded or die
	RESULT=$?
	if [ $RESULT -eq 0 ]; then
		echo "## SRA to FASTQ conversion succeeded for ${sra}"
	else
		echo "## SRA to FASTQ conversion failed for ${sra}, aborting!"
		exit 1
	fi

done

# return where you came from
cd -

# test for happy ending
if [ $? = 0 ]; then
echo "### all steps succeeded"
else
echo "### something ended wrong!, please check!"
fi