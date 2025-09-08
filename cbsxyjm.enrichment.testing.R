#very repetitive code to conduct enrichment testing using hotspot linkages 


library(gprofiler2)
library(ggplot2)

results = readRDS("/Users/noahalexander/s.cer.ensembl.112.RDS")


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrI")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrI"
hrange=c(1) #for the first/only hotspot on chrXVI
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")





#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrII")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrII"
hrange=c(1) #for the first/only hotspot on chrXVI
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrIV")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")


#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$Stress
dim(df)
df = subset(df, df$chr == "chrIV")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "Stress"
chrom.number = "chrIV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrIX")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrIX"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#


#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrVII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrVIII")
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

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrVIII"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(4)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrX"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#rest of chrx needs to be done

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXIII")
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


cross="3004"
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

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$S
dim(df)
df = subset(df, df$chr == "chrXV")
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


cross="3004"
experiment = "nacl"
subset = "S"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$Stress
dim(df)
df = subset(df, df$chr == "chrXV")
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


cross="3004"
experiment = "nacl"
subset = "Stress"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$Stress
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "Stress"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#

#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXVI")
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


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXVI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G2_M_mating"
chrom.number = "chrXVI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrXVI")
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


cross="3004"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrXVI"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G1_S
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
nn = names(table(df$bin))
cat(nn, sep=",")

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
nn = names(table(df$bin))
cat(nn, sep=",")
dim(df)

ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()


cross="3004"
experiment = "nacl"
subset = "G1_S"
chrom.number = "chrX"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
#

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
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


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange= hrange
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

#
