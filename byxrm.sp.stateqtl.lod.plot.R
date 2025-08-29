#BYxRM sp t0
df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/cell_cycle_assignment_LOD.RDS")

#df = df/(2*log(10))
df[1:6,1:10]

par(mar = c(5, 6, 4, 2) + 0.1)

plot(df[1,],                           
  type = "l",
  col = "light green",
  ylim = c(0, 50),
  ylab = "LOD",
  cex.lab = 2,
  cex.axis = 2,
  xaxt = "n",    
  xlab = "")

lines(df[2,],                             
  type = "l",
  col = "dark green")
lines(df[3,],                             
  type = "l",
  col = "light blue")
lines(df[4,],                             
  type = "l",
  col = "dark blue")
lines(df[5,],                             
  type = "l",
  col = "black")
lines(df[6,],                             
  type = "l",
  col = "dark red")



#this is to get the chr boundaries, the use of the esr vs pca-based state qtl files doesn't matter since the markers are shared 
df.t0 = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t0/state_pheno_LODs.RDS")
df.t0 = df.t0$LOD
df.t0 = as.data.frame(t(df.t0))
df.t0$observation <- 1:nrow(df.t0)
df.t0 = t(df.t0)
dim(df.t0)


dl = as.data.frame(t(df.t0))
dl$chr = word(rownames(dl),1,sep = "\\_")
head(dl)
table(dl$chr)


abline(a=NULL,b=NULL,h=4, col="red")
abline(a=NULL,b=NULL,h=NULL,v=c(0,cumsum(rle(dl$chr)$lengths)[1:15]),lty=2, col="black")



legend("topright",  box.lwd = 1,                      
  c("I", "II", "III", "IV", "V", "VI"),
  lty = 1,
  col = c("light green", "dark green", "light blue", "dark blue", "black", "dark red"))

axis(1, at = c(0,cumsum(rle(dl$chr)$lengths)[1:15]), labels = c(rle(dl$chr)$values), las=2)



################3 BYxRM sp t10 

df = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/cell_cycle_assignment_LOD.RDS")

#df = df/(2*log(10))
df[1:6,1:10]

par(mar = c(5, 6, 4, 2) + 0.1)

plot(df[1,],                           
  type = "l",
  col = "light green",
  ylim = c(0, 50),
  ylab = "LOD",
  cex.lab = 2,
  cex.axis = 2,
  xaxt = "n",    
  xlab = "")

lines(df[2,],                             
  type = "l",
  col = "dark green")
lines(df[3,],                             
  type = "l",
  col = "light blue")
lines(df[4,],                             
  type = "l",
  col = "dark blue")
lines(df[5,],                             
  type = "l",
  col = "black")
lines(df[6,],                             
  type = "l",
  col = "dark red")
lines(df[7,],                             
  type = "l",
  col = "dark red")



#this is to get the chr boundaries, the use of the esr vs pca-based state qtl files doesn't matter since the markers are shared 
df.t10 = readRDS("/Users/noahalexander/repeat_fine_mapping/combined/A/sp/t10/state_pheno_LODs.RDS")
df.t10 = df.t10$LOD
df.t10 = as.data.frame(t(df.t10))
df.t10$observation <- 1:nrow(df.t10)
df.t10 = t(df.t10)
dim(df.t10)


dl = as.data.frame(t(df.t10))
dl$chr = word(rownames(dl),1,sep = "\\_")
head(dl)
table(dl$chr)


abline(a=NULL,b=NULL,h=4, col="red")
abline(a=NULL,b=NULL,h=NULL,v=c(0,cumsum(rle(dl$chr)$lengths)[1:15]),lty=2, col="black")



legend("top",  box.lwd = 1,                      
  c("I", "II", "III", "IV", "V", "VI", "VII"),
  lty = 1,
  col = c("light green", "dark green", "light blue", "dark blue", "black", "dark red", "purple"))

axis(1, at = c(0,cumsum(rle(dl$chr)$lengths)[1:15]), labels = c(rle(dl$chr)$values), las=2)
