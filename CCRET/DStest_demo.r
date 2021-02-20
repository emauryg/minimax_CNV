##--------------------------------------------------
## Convert plink data format to CCRET data format
##...................................................
## Please see "readme.CCRET_data_conversion.txt" for
## detailed descriptions.
##--------------------------------------------------
source('fun.plink2ccret.r')
print(date()) #N=2000
plink2ccret(tag='mycnv1',  gset='mygeneset.bed', workdir='./Data/')



##------------------------------------
## Read in data (in ccret format)
##------------------------------------
tag    = "mycnv1"    
yy     = unlist(read.table(paste(tag, "_yy.txt",sep=''), header=T))  
ds.mat =        read.table(paste(tag, "_ds.txt",sep=''), header=T)
ln.mat =        read.table(paste(tag, "_ln.txt",sep=''), header=T)
gi.mat =        read.table(paste(tag, "_gi.txt",sep=''), header=T)

nsubj  = nrow(ds.mat); nsubj ## 2000
ncnvr  = ncol(ds.mat); ncnvr ## 1320    
dim(gi.mat)   #2000, 668

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
actual.freq = colMeans(ds.mat!=2)
key.keep    = (actual.freq>0)
ds.mat      = ds.mat[,key.keep]
ln.mat      = ln.mat[,key.keep]

nloci       = ncol(ds.mat)
nsubj       = nrow(ds.mat)
cat("# of loci=",nloci,"| # of subj=",nsubj,"\n")

##-----------------------------
##  calculating kernel matrix
##-----------------------------
dsnew    = ds.to.dsnew(ds.mat)
Kmat.ds  = tcrossprod(dsnew)
##-------------------------------------------------------------------
## obtain p-value of CCRET test for dosage effects, NOT adjusting for length and GI
##-------------------------------------------------------------------
pv.nocov = vctest.btqt.Gmain.fun(y=yy, geno=ds.mat, x.adj=NA, trait.type="binomial", SSS=Kmat.ds)

##-------------------------------------------------------------------
## obtain p-value of CCRET test for dosage effects, adjusting for length and GI
##-------------------------------------------------------------------
gi.indv             = rowSums(gi.mat!=0)
Ncnv.indv           = rowSums(ds.mat!=2)
raw.lavg.ln         = apply(ln.mat, 1, function(x){mean(x[x>0])}) ## nsubj
lavg.ln.overallmean = raw.lavg.ln; lavg.ln.overallmean[is.na(raw.lavg.ln)] =mean(ln.mat[ln.mat>0])
pv.withcov          = vctest.btqt.Gmain.fun(y=yy, geno=ds.mat, x.adj=cbind(gi.indv, Ncnv.indv), trait.type="binomial", SSS=Kmat.ds)

##-------------------------------------------------------------------------------
## output of DS analysis:
##        CCRET will output two p-values. Note that the first p-value is obtained
##        using Davies method and the second p-value is obtained using Liu et al.
##        method. (Please see more details in 'CompQuadForm' package).
##-------------------------------------------------------------------------------
cat("CCRET p-value of DS effect NOT adjusting for length and gene intersection=",pv.nocov,"\n")
################CCRET p-value of DS effect NOT adjusting for length and gene intersection= 9.563278e-06 2.192841e-06 

cat("CCRET p-value of DS effect adjusting for length and gene intersection=",pv.withcov,"\n")
#################CCRET p-value of DS effect adjusting for length and gene intersection= 0 5.712077e-07 



