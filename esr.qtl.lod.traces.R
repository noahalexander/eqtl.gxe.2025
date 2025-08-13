##################  BYxRM nacl experiment

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t0/state_pheno_LODs.RDS")
df = df$LOD
df = as.data.frame(t(df))
df$observation <- 1:nrow(df)
df = t(df)
dim(df)




plot(df[4,],                           
  df[1,],
  type = "l",
  ylim = c(0, 75),
  col = "light green",
  ylab = "LOD",
  cex.lab = 2,
  cex.axis = 2,
  xaxt = "n",    
  xlab = ""      
)
lines(df[4,],                             
  df[2,],
  type = "l",
  col = "dark green")
lines(df[4,],                             
  df[3,],
  type = "l",
  col = "light blue"
)


df.t30 = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/nacl/t30/state_pheno_LODs.RDS")
df.t30 = df.t30$LOD
df.t30 = as.data.frame(t(df.t30))
df.t30$observation <- 1:nrow(df.t30)
df.t30 = t(df.t30)
dim(df.t30)


lines(df.t30[4,],                             
  df.t30[1,],
  type = "l",
  col = "orange")
lines(df.t30[4,],                             
  df.t30[2,],
  type = "l",
  col = "red")
lines(df.t30[4,],                             
  df.t30[3,],
  type = "l",
  col = "yellow")

#the markers are the same in the t0 and t30 datasets
dl = as.data.frame(t(df.t30))
dl$chr = word(rownames(dl),1,sep = "\\_")
head(dl)
table(dl$chr)



abline(h=4, col="red")
abline(a=NULL,b=NULL,h=NULL,v=c(0,cumsum(rle(dl$chr)$lengths)[1:15]),lty=2, col="black")
#abline(a=NULL,b=NULL,h=NULL,v=test.vec,lty=2, col="red")


legend("topleft",                           
  c("t0 RP", "t0 iESR", "t0 RiBi", "t30 RP", "t30 iESR", "t30 RiBi"),
  lty = 1,
  col = c("light green", "dark green", "light blue", "orange", "red", "yellow"))

axis(1, at = c(0,cumsum(rle(dl$chr)$lengths)[1:15]), labels = c(rle(dl$chr)$values), las=2)

