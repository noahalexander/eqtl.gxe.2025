library(gprofiler2)
library(ggplot2)

results = readRDS("/Users/noahalexander/s.cer.ensembl.112.RDS")
setwd("/Users/noahalexander/Documents/enrichment.output.20241030/")


enrich.cust <- function(cross,experiment, timepoint, subset, chrom.number, hrange, thresh = 0.05, correction_meth = "g_SCS") {
  #function to generate gprofiler-based GO enrichment results for transcripts linking to hotspots 
  #elist.out[[1]]$result is for pos betas 
  #elist.out[[2]]$result is for neg betas
  #sometimes there will only be one or none depending on whether the input sets of transcripts have enriched terms 
  #example inputs:
  #cross="A"
  #experiment = "sp"
  #subset = "IV"
  #chrom.number = "chrXV"
  #hrange=c(1:2) for first two hotspots on chrXV
  #thresh = 0.05
  #correction_meth = "g_SCS"
  
  #subset biomart output to the chromosome where the hotspot is using existing 'results' object/the output of biomart getBM
  dd = subset(results, results$chromosome_name == substr(chrom.number, 4, nchar(chrom.number)))
  #prepare custom background of all genes used for eqtl mapping (segData$Yr), excluding the genes on the same chromosome as the hotspot
  segdata= readRDS(paste0("/Users/noahalexander/repeat_fine_mapping/combined/",cross,"/",experiment,"/",timepoint,"/segData.RDS"))
  background = colnames(segdata$Yr)
  length(background)
  background = subset(background, !(background %in% dd$ensembl_gene_id))
  length(background)
  
  df = readRDS(paste0("/Users/noahalexander/repeat_fine_mapping/combined/",cross,"/",experiment,"/",timepoint,"/hotspot_peaks.RDS"))
  df = df[[subset]]
  dim(df)
  df = subset(df, df$chr == chrom.number)
  dim(df)
  df = subset(df, df$in.hotspot == "TRUE")
  dim(df)
  unique(df$bin)
  table(df$bin)
  
  #subset for a range of bins (where 1 is the bin closest to the start of the chromosome)
  df = subset(df, df$bin %in% names(table(df$bin))[hrange])
  
  
  gostres.total <- gost(query = df$transcript, 
    organism = "scerevisiae", ordered_query = FALSE, 
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
    measure_underrepresentation = FALSE, evcodes = FALSE, 
    user_threshold = thresh, correction_method = correction_meth, 
    domain_scope = "annotated", custom_bg = background, 
    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
  head(gostres.total$result)
  gostres.total$result$parents <-NULL
  
  if (!is.null(gostres.total)) {
    gostres.total$result <- gostres.total$result[order(gostres.total$result$p_value),]
  }

  
  
  filename = paste0(cross,".", experiment,".", timepoint,".", subset,".", chrom.number,".", paste0(unique(df$bin), collapse = "_"), ".",thresh,".", correction_meth,".", "gostres.total.csv") 
  write.table(as.data.frame(gostres.total$result), file=filename, quote=F, sep=",", row.names=F, col.names=T)
  
  
  df.pos = subset(df, df$Beta >0)
  df.neg = subset(df, df$Beta <0)
  
  gostres.pos <- gost(query = df.pos$transcript, 
    organism = "scerevisiae", ordered_query = FALSE, 
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
    measure_underrepresentation = FALSE, evcodes = FALSE, 
    user_threshold = thresh, correction_method = correction_meth, 
    domain_scope = "annotated", custom_bg = colnames(segdata$Yr), 
    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
  head(gostres.pos$result)
  gostres.pos$result$parents <-NULL
  
  if (!is.null(gostres.pos)) {
  gostres.pos$result <- gostres.pos$result[order(gostres.pos$result$p_value),]
  }
  
  filename = paste0(cross,".", experiment,".", timepoint,".", subset,".", chrom.number,".", paste0(unique(df.pos$bin), collapse = "_"), ".",thresh,".", correction_meth,".", "gostres.pos.csv") 
  write.table(as.data.frame(gostres.pos$result), file=filename, quote=F, sep=",", row.names=F, col.names=T)
  
  
  gostres.neg <- gost(query = df.neg$transcript, 
    organism = "scerevisiae", ordered_query = FALSE, 
    multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
    measure_underrepresentation = FALSE, evcodes = FALSE, 
    user_threshold = thresh, correction_method = correction_meth, 
    domain_scope = "annotated", custom_bg = colnames(segdata$Yr), 
    numeric_ns = "", sources = NULL, as_short_link = FALSE, highlight = TRUE)
  head(gostres.neg$result)
  gostres.neg$result$parents <-NULL


  if (!is.null(gostres.neg)) {
    gostres.neg$result <- gostres.neg$result[order(gostres.neg$result$p_value),]
  }
  
  filename = paste0(cross,".", experiment,".", timepoint,".", subset,".", chrom.number,".", paste0(unique(df.neg$bin), collapse = "_"), ".",thresh,".", correction_meth,".", "gostres.neg.csv")
  write.table(as.data.frame(gostres.neg$result), file=filename, quote=F,sep=",", row.names=F, col.names=T)
  
  
  filename = paste0(cross,".", experiment,".", timepoint,".", subset,".", chrom.number,".", paste0(unique(df$bin), collapse = "_"), ".",thresh,".", correction_meth,".", "beta_se.png")
  ggplot(df, aes(x=df$Beta/df$SE)) + geom_histogram()
  ggsave(filename, device = "png")
  
  elist = list(c(gostres.pos, gostres.neg))
  
  return(elist)
}
