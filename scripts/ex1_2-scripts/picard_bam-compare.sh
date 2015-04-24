#!/bin/bash
# picard_bam-compare.sh

# compair 'all' mappings pairwise

results=trimming-effect_SRR1039509/picard_bam-compare
mkdir -p ${results}

# parameters
orimap="tophat2.0.13-SRR1039509_results/SRR1039509_all-tophat.bam"
trimmedmap="tophat2.0.13-SRR1039509-trimmed_results/SRR1039509_all-tophat.bam"

java -Xmx4g -jar $PICARD/picard.jar CompareSAMs \
	${orimap} \
	${trimmedmap} \
	>${results}/picard_compare_${met}.txt \
	2>${results}/CompareSAMs_${met}.log
