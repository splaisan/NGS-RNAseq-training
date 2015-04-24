#!/usr/bin/RScript
# compare-htseq_counts.R
# copy this code to RStudio and adapt file locations to match yours

basedir <- "/media/bits/RNASeq_DATA/htseq_counts"
setwd(basedir)

# Take htseq-count for 'all' and 'gtf' and make scatter-plot
# use sample SRR1039509 for plot
SRR1039509_all <- read.delim("SRR1039509_all_counts.txt", header=F)
SRR1039509_gtf <- read.delim("SRR1039509_gtf_counts.txt", header=F)

data <- merge(SRR1039509_all, SRR1039509_gtf, by='V1')
colnames(data) <- c("ID", "SRR1039509_all", "SRR1039509_gtf")
head(data, 10)

# ID SRR1039509_all SRR1039509_gtf
# 1  __alignment_not_unique        5535220        2553851
# 2             __ambiguous         567580         691622
# 3            __no_feature        1871251        1212702
# 4           __not_aligned              0              0
# 5         __too_low_aQual              0              0
# 6         ENSG00000000003            434            440
# 7         ENSG00000000005              0              0
# 8         ENSG00000000419            488            506
# 9         ENSG00000000457            226            234
# 10        ENSG00000000460             52             56

# extract summary lines
data.summary <- data[grep("^ENS", data$ID, perl=TRUE, invert=TRUE), ]
rownames(data.summary) <- data.summary$ID
data.summary <- data.all.summary[,-1]

data.summary

# ID SRR1039509_all SRR1039509_gtf
# __alignment_not_unique __alignment_not_unique        5535220        2553851
# __ambiguous                       __ambiguous         567580         691622
# __no_feature                     __no_feature        1871251        1212702
# __not_aligned                   __not_aligned              0              0
# __too_low_aQual               __too_low_aQual              0              0

# extract gene data lines
data <- data[grep("^ENS", data$ID, perl=TRUE, invert=FALSE), ]
head(data)

# ID SRR1039509_all SRR1039509_gtf
# 6  ENSG00000000003            434            440
# 7  ENSG00000000005              0              0
# 8  ENSG00000000419            488            506
# 9  ENSG00000000457            226            234
# 10 ENSG00000000460             52             56
# 11 ENSG00000000938              0              0

# compute correlation coefficient between both datasets
cor(data$SRR1039509_all, data$SRR1039509_gtf, 
    use="complete.obs", method="kendall")

# [1] 0.9230123

# correlate with more control
cor.test(data$SRR1039509_all, data$SRR1039509_gtf,
         alternative = "two.sided",
         method = "kendall",
         exact = NULL, 
         conf.level = 0.95)

# Kendall's rank correlation tau
# 
# data:  data$SRR1039509_all and data$SRR1039509_gtf
# z = 266.8679, p-value < 2.2e-16
# alternative hypothesis: true tau is not equal to 0
# sample estimates:
#   tau 
# 0.9230123 

# plot scatterplot
# png(file="SRR1039509-all_vs_gtf-counts.png")

plot(log2(data$SRR1039509_all) ~ log2(data$SRR1039509_gtf),
     xlab="log2-GTF counts - SRR1039509",
     ylab="log2-ALL counts - SRR1039509",
     pch=20,
     cex=0.5
  )

# add lines
abline(0, 1, col="red", lty=1, lwd=2)
abline(h=log2(10), col="grey1", lty=2)
abline(v=log2(10), col="grey1", lty=2)
abline(h=log2(20), col="grey2", lty=2)
abline(v=log2(20), col="grey2", lty=2)
# dev.off()
