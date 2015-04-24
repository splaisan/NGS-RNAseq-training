# define working folder
workdir <- "/media/bits/RNASeq_DATA"
setwd(workdir)

# get information from the ebi
PID <- "SRP033351"
ena.url <- paste("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
                 PID,
                 "&result=read_run",
                 "&fields=run_accession,library_name,",
                 "read_count,fastq_ftp,fastq_aspera,",
                 "fastq_galaxy,sra_ftp,sra_aspera,sra_galaxy,",
                 "&download=text",
                 sep="")
metadata <- read.table(url(ena.url), header=TRUE, sep="\t")

# get first line and transpose
t(metadata[1,])

# construct metadata table from the online information

library("stringr")

PID <- "SRP033351"
ena.url <- paste("http://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=",
                 PID,
                 "&result=read_run",
                 "&fields=run_accession,read_count",
                 sep="")

# assemble into a data.frame
metadata <- as.data.frame(read.table(url(ena.url), header=TRUE, sep="\t"))

# add missing information about cells & treatment
samples=c("N61311_untreated",
          "N61311_Dex",
          "N61311_Alb",
          "N61311_Alb_Dex",
          "N052611_untreated",
          "N052611_Dex",
          "N052611_Alb",
          "N052611_Alb_Dex",
          "N080611_untreated",
          "N080611_Dex",
          "N080611_Alb",
          "N080611_Alb_Dex",
          "N061011_untreated",
          "N061011_Dex",
          "N061011_Alb",
          "N061011_Alb_Dex")
metadata$samples <- samples

# sample the first line
t(metadata[1,])

# add new columns
metadata$cells <- str_extract(metadata$sample, "\\b[A-Za-z0-9]+")
metadata$treatment <- mapply(sub,paste(metadata$cells,"_",sep=""), "", metadata$sample)

# view table
metadata

# save to file for reuse
write.table(metadata, file="GSE52778_metadata.txt", row.names=F, sep="\t", quote=FALSE)
