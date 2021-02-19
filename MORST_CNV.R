
## R code to run MORST_CNV, based on  great part from https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1831926

## Written by Eduardo Maury (eduardo_maury AT hms DOT harvard DOT edu)

library("MASS")
library("CompQuadForm")
library("matrixStats")
library("Matrix")
library("Rcpp")
library("RcppEigen")
sourceCpp("Saddle.cpp")
source("Association_tests.R")


## Preprocessing as recommended in the CCRET package
##------------------------------------
## Read in data (in ccret format)
##------------------------------------
tag    = "data/mycnv2" 
yy     = unlist(read.table(paste(tag, "_yy.txt",sep=''), header=T))
ds.mat =        read.table(paste(tag, "_ds.txt",sep=''), header=T)
ln.mat =        read.table(paste(tag, "_ln.txt",sep=''), header=T)
gi.mat =        read.table(paste(tag, "_gi.txt",sep=''), header=T)


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
n.cnv = rowSums(ds.mat!=2)
raw.lavg.ln = apply(ln.mat, 1, function(x){mean(x[x>0])}) ## nsubj
lavg.ln.overallmean = raw.lavg.ln
lavg.ln.overallmean[is.na(raw.lavg.ln)] =mean(ln.mat[ln.mat>0])

## Create fixed-design matrix
## Note that it gets collapsed into a n x 2 matrix. 
Z=cbind(n.cnv,lavg.ln.overallmean)

## Calculate null model 
Y = yy

obj <- MORST_NULL_Model(Y, Z)


X = as.matrix(gi.mat)

morst_res = MORST(G=X, obj=obj, weights= rep(1/ngene,ngene), siglevel = 0.05)

print(paste("MORST P-Value:",morst_res))
