library(GenomicRanges)
library(IRanges)
library(stringr)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)
gr.hotspot.bins.3004.nacl.t30 = with(dd, GRanges(chr, IRanges(start=start, end = end), freq=freq))


write.table(dd, file="3004.t30.combined.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$Stress
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)
gr.hotspot.bins.3004.nacl.t30 = with(dd, GRanges(chr, IRanges(start=start, end = end), freq=freq))


write.table(dd, file="3004.t30.stress.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t30/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)
gr.hotspot.bins.3004.nacl.t30 = with(dd, GRanges(chr, IRanges(start=start, end = end), freq=freq))


write.table(dd, file="3004.t30.s.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)





df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0/hotspot_peaks.RDS")
df = df$G2_M_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3004.t0.g2m.mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)





df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.combined.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$I
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.I.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$II
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.II.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$III
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.III.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$IV
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.IV.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$V
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.V.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/hotspot_peaks.RDS")
df = df$VI
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t0.VI.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)




















df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.combined.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)



df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$I
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.I.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)




df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$II
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.II.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$III
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.III.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$IV
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.IV.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)



df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$V
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.v.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$VI
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.VI.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/hotspot_peaks.RDS")
df = df$VII
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.sp.t10.VII.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)



##################salt 

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.combined.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)



df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$G1
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.g1.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$G1_S
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.g1.s.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$G2_M_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.G2_M_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$G2_M_no_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.G2_M_no_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)



df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$M_G1_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.M_G1_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t0/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t0.s.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)





df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$combined
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.combined.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$G1
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.g1.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$G1_S
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.g1.s.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$G2_M_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.G2_M_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$G2_M_no_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.G2_M_no_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$M_G1_no_mating
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.M_G1_no_mating.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)


df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl//t30/hotspot_peaks.RDS")
df = df$S
df = subset(df, df$in.hotspot == TRUE)
dd = data.frame(chr = df$chrom, start = df$CI.l, end=df$CI.r, LOD=df$LOD, bin=df$bin)
unique_bins = unique(dd$bin)
bin_table = table(dd$bin)
bin_table = as.data.frame(bin_table)
dd = data.frame(chr=word(unique_bins,1,sep = ":"), start=as.numeric(word(word(unique_bins,1,sep = "-"), 2, sep = ":")), end = as.numeric(word(unique_bins,2,sep = "\\-")), freq=bin_table$Freq)

write.table(dd, file="3051.nacl.t30.S.hotspots.bed", quote=F, sep="\t", row.names=F, col.names=T)
