##--------------------------------------------------
## Convert plink data format to CCRET data format
##...................................................
## Please see "readme.CCRET_data_conversion.txt" for
## detailed descriptions.
##--------------------------------------------------
#setwd("C:/Users/emaur/Dropbox (MIT)/PhD Classes/STAT364 Scalable Statistical Inference for Big Data and Applications/presentation_Lin/CCRET_code/")
source('fun.plink2ccret.r')
print(date()) #N=2000
plink2ccret(tag='mycnv2',  gset='mygeneset.bed', workdir='./Data/')

##------------------------------------
## Read in data (in ccret format)
##------------------------------------
tag    = "mycnv2" 
yy     = unlist(read.table(paste(tag, "_yy.txt",sep=''), header=T))
ds.mat =        read.table(paste(tag, "_ds.txt",sep=''), header=T)
ln.mat =        read.table(paste(tag, "_ln.txt",sep=''), header=T)
gi.mat =        read.table(paste(tag, "_gi.txt",sep=''), header=T)

##------------------------------------------------
## source the functions needed for running CCRET
##------------------------------------------------
getwd()
setwd("../") ## this change the working directory to one level abouve /Data/
getwd()
source("corefun.r")

##---------------------------------
## only keep loci with >0 CNV event 
##---------------------------------
actual.freq = colMeans(gi.mat!=0)
key.keep    = (actual.freq>0)
gi.mat      = gi.mat[,key.keep]
ngene       = ncol(gi.mat)
nsubj       = nrow(gi.mat)
##-----------------------
## calculating avg len and CNV count (n.cnv)
##-----------------------
n.cnv               = rowSums(ds.mat!=2)
raw.lavg.ln         = apply(ln.mat, 1, function(x){mean(x[x>0])}) ## nsubj
lavg.ln.overallmean = raw.lavg.ln; lavg.ln.overallmean[is.na(raw.lavg.ln)] =mean(ln.mat[ln.mat>0])
##----------------------------------------
## calculating kernel matrix based on GI
## obtain CCRET p-value
##----------------------------------------
GC1mat = getK.gi.mat.fun(gi.mat, kernel=1, WT=rep(1/ngene, ngene)) ## linear kernel
pv     = vctest.btqt.Gmain.fun(y=yy, geno=gi.mat, x.adj=cbind(n.cnv,lavg.ln.overallmean), trait.type="binomial", SSS=GC1mat)

##-------------------------------------------------------------------------------
## output of GI analysis
##        CCRET will output two p-values. Note that the first p-value is obtained
##        using Davies method and the second p-value is obtained using Liu et al.
##        method. (Please see more details in 'CompQuadForm' package).
##-------------------------------------------------------------------------------
cat("CCRET p-value of GI effect adjusting for length and dosage =",pv,"\n")
######### CCRET p-value of GI effect adjusting for length and dosage = 0 6.109567e-09 


##----------------------------
##  fixed effect analysis 
##----------------------------
##  ngi = rowSums(gi.mat!=0); out=summary(glm(yy~ngi+n.cnv+lavg.ln.overallmean)); print(out$coef)

## ##---
## ##                         Estimate   Std. Error   t value     Pr(>|t|)
## ## (Intercept)         4.328319e-01 2.438153e-02 17.752449 1.389745e-65
## ## ngi                 1.063568e-01 3.353941e-02  3.171100 1.541684e-03
## ## n.cnv               1.079437e-02 6.931860e-03  1.557211 1.195789e-01
## ## lavg.ln.overallmean 2.571539e-07 1.220839e-07  2.106371 3.529671e-02

