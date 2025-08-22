
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/3004/nacl/t0//segData.RDS")

Gsub = df$Gsub
#start with column names of gsub and strip out the sequence position
gcols = colnames(Gsub)
index = data.frame(chr=word(gcols,1,sep = "\\_"), start= word(gcols,2,sep = "\\_"))
#make a granges using the positions of each marker (not ranges but positions). 
gr.index= with(index, GRanges(chr, IRanges(start=start)))

#section 5"Finding the nearest genomic position in GRanges objects" in https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
#for example, this granges object is 'x' and the index of markers from Gsub is 'subject'
nearest.markers = nearest(gr.ie.3004.nacl.t0, gr.index)

#make new set of start and end positions using nearest/the index generated 
####the naming of these functions is odd, follow and precede give rightmost and leftmost adjacent markers it seems... https://bioconductor.org/packages/devel/bioc/vignettes/GenomicRanges/inst/doc/GenomicRangesIntroduction.html
# vector of indices that follow the interval
right.markers = precede(gr.ie.3004.nacl.t0, gr.index)

# vector of indices that preced the interval
left.markers = follow(gr.ie.3004.nacl.t0, gr.index)

gr.ie.3004.nacl.t0$nearest = nearest.markers
gr.ie.3004.nacl.t0$right.markers = right.markers
gr.ie.3004.nacl.t0$left.markers = left.markers
  
#figure out why NAs happen/where
#dd= gr.ie.3004.nacl.t0
#subset(dd, is.na(dd$right.markers == T))

#isub.right = index[right.markers,]
#isub.left = index[left.markers,]





#############---------------#################-------------------###################################3004 cross NaCl timecourse 
#read in segdata and make granges object using gsub. this can be used for all time points in this experiment   

sdf = readRDS("/Users/noahalexander/eqtl.gxe/results/combined/3004.nacl.t0.combined.out.20230921/segData.RDS")
Gsub = sdf$Gsub
#start with column names of gsub and strip out the sequence position
gcols = colnames(Gsub)
index = data.frame(chr=word(gcols,1,sep = "\\_"), start= word(gcols,2,sep = "\\_"))
#make a granges using the positions of each marker (not ranges but positions). 
gr.index= with(index, GRanges(chr, IRanges(start=start)))


####---------------------------------3004
#3004 salt stimulation t0

#start with state qtl peaks
df = readRDS("/Users/noahalexander/eqtl.gxe/results/combined/3004.nacl.t0.combined.out.20230921/state_pheno_LOD_peaks.RDS")
#subset for rp assignments 
dfrp.3004.nacl.t0 = df[[1]]
#subset for iesr assignments 
#dfie.3004.nacl.t0 = df[[2]]
#subset for ribi assignments 
#dfrb.3004.nacl.t0 = df[[3]]

##make data frame for each state and make a granges 
drp.3004.nacl.t0 = data.frame(chr=word(dfrp.3004.nacl.t0$fscan.markers,1,sep = "\\_"), 
  start= sapply(strsplit(dfrp.3004.nacl.t0$CI.l, "_"), function(x) x[2]), 
  end= sapply(strsplit(dfrp.3004.nacl.t0$CI.r, "_"), function(x) x[2] ),
  p=dfrp.3004.nacl.t0$p,
  q=dfrp.3004.nacl.t0$q,
  LOD=dfrp.3004.nacl.t0$LOD,
  fscan.markers= dfrp.3004.nacl.t0$fscan.markers
)
drp.3004.nacl.t0$start = as.numeric(drp.3004.nacl.t0$start)
#rangesB$start= as.numeric(rangesB$start)
drp.3004.nacl.t0$end = as.numeric(drp.3004.nacl.t0$end)

#make granges object
gr.rp.3004.nacl.t0= with(drp.3004.nacl.t0, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
gr.rp.3004.nacl.t0= gr.rp.3004.nacl.t0[order(gr.rp.3004.nacl.t0$LOD, decreasing = T),]

#look for markers from index granges that are to the left or right of the input interval (state qtl CI)
#the output is an index that can be used to search the index granges for the adjacent markers to widen the interval 
#these are indices, not bp coordinates but they can be converted using index granges 

nearest.markers = nearest(gr.rp.3004.nacl.t0, gr.index)
right.markers = precede(gr.rp.3004.nacl.t0, gr.index)

# vector of indices that precede the interval in spite of the function name 
left.markers = follow(gr.rp.3004.nacl.t0, gr.index)
#add positions of markers on either side of returned interval 
gr.rp.3004.nacl.t0$left.markers = left.markers
gr.rp.3004.nacl.t0$nearest = nearest.markers
gr.rp.3004.nacl.t0$right.markers = right.markers


#pull out the bp coordinates using the indices 
isub.left = index[left.markers,]
left.coorinate = as.numeric(isub.left$start)

isub.nearest = index[nearest.markers,]
nearest.coordinate = as.numeric(isub.nearest$start)

isub.right = index[right.markers,]
right.coorinate = as.numeric(isub.right$start)

#add bp coordinates on either side of Josh's returned interval 
gr.rp.3004.nacl.t0$left.coorinate = left.coorinate
gr.rp.3004.nacl.t0$nearest.coordinate = nearest.coordinate
gr.rp.3004.nacl.t0$right.coorinate = right.coorinate

#if not containing an na, the difference should reflect the size of the new CI
gr.rp.3004.nacl.t0$new_CI_length = right.coorinate - left.coorinate

#deal with NAs and remake ranges to reflect the new interval
df = as.data.frame(gr.rp.3004.nacl.t0)
for(i in 1:nrow(df)) {
  row <- df[i,]
  #print(row[2])
  if (is.na(row[12]) == T) {
    #set left.marker to the original start of the interval
    df[i,12] = df[i,2]
    #print(row[9])
  }
  if (is.na(row[14]) == T) {
    #set right.marker to the original end of the interval
    df[i,14] = df[i,3]
    #print(row[9])
  }
}

#remake granges adding the new/checked broader CIs
#gr.rp.3004.nacl.t0= with(drp.3004.nacl.t0, GRanges(chr, IRanges(start=start, end=end), LOD=LOD, p=p, q=q))
#using column names for new data frame, make granges with wider range
gr.rp.3004.nacl.t0 = with(df, GRanges(seqnames = seqnames, IRanges(start=left.coorinate, end=right.coorinate), LOD=LOD, p=p, q=q, left.markers = left.markers, nearest.markers=nearest, right.markers=right.markers, original.start=start, original.end=end))

as.data.frame(gr.rp.3004.nacl.t0)


#######testing using fxn 


dd = expand_interval(gr.state = gr.ie.3004.nacl.t0, gr.experiment.index = gr.index)
