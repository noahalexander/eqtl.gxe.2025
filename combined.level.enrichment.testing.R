#very repetitive code to conduct enrichment testing using hotspot linkages 


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrI"
hrange=c(1) #for the first/only hotspot on chrXVI
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrII"
hrange=c(2,3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
######--------------chrIII
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIII"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrIX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrIX"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrVII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrVII"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrVII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrVII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrVIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrVIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrVIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrX"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrX"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrX"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrX")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrX"
hrange=c(2,3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXII"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXII"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXII"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##
##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXII"
hrange=c(3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXII"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")

##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXIII"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIII")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIII"
hrange=c(2,3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2,3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(2,3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)


cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXIV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(3)
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXIV"
hrange=c(3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,4) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(1,4) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange=c(1,2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1,3) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(1,3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(2) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange=c(2) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(3) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXV"
hrange=c(3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXV")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(3) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXV"
hrange=c(3) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXVI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="3004"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXVI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXVI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t0"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "sp"
subset = "combined"
chrom.number = "chrXVI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t10"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##

##
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/hotspot_peaks.RDS")
df = df$combined
dim(df)
df = subset(df, df$chr == "chrXVI")
dim(df)
df = subset(df, df$in.hotspot == "TRUE")
dim(df)
names(table(df$bin))

hrange=c(1) 
names(table(df$bin))[hrange]
df = subset(df, df$bin %in% names(table(df$bin))[hrange])
names(table(df$bin))
dim(df)

cross="A"
experiment = "nacl"
subset = "combined"
chrom.number = "chrXVI"
hrange=c(1) 
thresh = 0.05
correction_meth = "g_SCS"
timepoint = "t30"
outlist = enrich.cust(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS")
##
