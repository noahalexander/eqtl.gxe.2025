#very repetitive code to conduct enrichment testing using hotspot linkages 

library(gprofiler2)
library(ggplot2)

results = readRDS("/Users/noahalexander/s.cer.ensembl.112.RDS")

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#for the adjacent hotspot
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrX"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrX"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(3,4)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$IV
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "IV"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "I"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "III"
chrom.number = "chrXVI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "sp"
subset = "II"
chrom.number = "chrXVI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#################-------------------------------------------------------------------------------NaCl











#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)
ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()

cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)
ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()

cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)
ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()

cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_no_mating"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrVIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrVIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(4)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_no_mating"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$M_G1_no_mating
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "M_G1_no_mating"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_no_mating"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$M_G1_no_mating
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "M_G1_no_mating"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$M_G1_no_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "M_G1_no_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_no_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_no_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "S"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#
#
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$G1
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

#for the leftmost 
hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="A"
experiment = "nacl"
subset = "G1"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
#

