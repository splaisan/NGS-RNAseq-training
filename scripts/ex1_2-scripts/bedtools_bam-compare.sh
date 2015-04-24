#!/bin/bash
# bedtools_bam-compare.sh

# compare the results of 'all' tophat mapping
# between unprocessed and trimmed SRR1039509 reads

# parameters
orimap="tophat2.0.13-SRR1039509_results/SRR1039509_all-tophat.bam"
trimmedmap="tophat2.0.13-SRR1039509-trimmed_results/SRR1039509_all-tophat.bam"

results=trimming-effect_SRR1039509/bedtools_bam-compare
mkdir -p ${results}

# extract the coverage from each bam in bedgraph format file using bedtools genomecov
bedtools genomecov -bg \
    -ibam ${orimap} \
    | bgzip -c > ${results}/orimap.bedgraph.gz

bedtools genomecov -bg \
    -ibam ${trimmedmap} \
    | bgzip -c > ${results}/trimmedmap.bedgraph.gz

# process both coverage bedgraph files using bedtools unionbedg
hg19idx=$BOWTIE2_INDEXES/hg19.fa.fai

bedtools unionbedg -header \
    -i ${results}/orimap.bedgraph.gz ${results}/trimmedmap.bedgraph.gz \
    -names 'raw-reads' 'trimmed-reads' \
    -g ${hg19idx} \
    -empty \
    | bgzip -c > ${results}/bedtools-unionbedg.bg.gz
