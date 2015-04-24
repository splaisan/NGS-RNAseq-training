#!/usr/bin/RScript
# htseq-combine_all.R
 
# Take all htseq-count results and melt them in to one big dataframe
## do this for either tophat_all or tophat_gtf mappings.
basedir <- "/media/bits/RNASeq_DATA"
setwd(basedir)
 
rawfiles <- paste(basedir, "htseq_counts", sep="/")
tophat.all <- list.files(path = rawfiles, 
	pattern = "_all_counts.txt",
	all.files = TRUE, 
	recursive = FALSE, 
	ignore.case = FALSE, 
	include.dirs = FALSE)

# we choose the 'all' series
myfiles <- tophat.all
DT <- list()
 
# read each file as array element of DT and rename the last 2 cols
for (i in 1:length(myfiles) ) {
  infile = paste("htseq_counts", myfiles[i], sep = "/")
	DT[[myfiles[i]]] <- read.table(infile, header = F)
	cnts <- gsub("(.*)_all_mdup_counts.txt", "\\1", myfiles[i])
	colnames(DT[[myfiles[i]]]) <- c("ID", cnts)
}
 
# merge all elements based on first ID columns
data <- DT[[myfiles[1]]]
for (i in 2:length(myfiles)) {
	y <- DT[[myfiles[i]]]
	z <- merge(data, y, by = c("ID"))
	data <- z
}
 
# remove top two rows
data.all.summary <- data[grep("^ENS", data$ID, perl=TRUE, invert=TRUE), ]
rownames(data.all.summary) <- data.all.summary$ID
data.all.summary <- data.all.summary[,-1]
 
data.all <- data[grep("^ENS", data$ID, perl=TRUE, invert=FALSE), ]
 
# final merged table
head(data.all, 3)
 
# summary counts by type
t(data.all.summary)
 
write.table(data.all,
	file = "htseq-counts_all.csv",
	append = FALSE,
	quote = FALSE,
	sep = " ",
	eol = "\n",
	na = "NA",
	dec = ".",
	row.names = FALSE,
	col.names = TRUE,
	qmethod = c("escape", "double"))
 
write.table(data.all.summary,
	file = "htseq-counts_all-summary.csv",
	append = FALSE,
	quote = FALSE,
	sep = " ",
	eol = "\n",
	na = "NA",
	dec = ".",
	row.names = FALSE,
	col.names = TRUE,
	qmethod = c("escape", "double"))
 
# cleanup intermediate objects
rm(y, z, i, data)
