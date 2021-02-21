
## R code to run simulations comparing CCRET and MORST approach to CNV geneset association

## Written by Eduardo Maury (eduardo_maury AT hms DOT harvard DOT edu)

setwd("C:/Users/emaur/Dropbox (MIT)/PhD Classes/STAT364 Scalable Statistical Inference for Big Data and Applications/presentation_Lin/minimax_CNV")

library("MASS")
library("CompQuadForm")
library("matrixStats")
library("Matrix")
library("Rcpp")
library("RcppEigen")
sourceCpp("Saddle.cpp")
source("Association_tests.R")
source("CCRET/corefun.r")
source('CCRET/fun.plink2ccret.r')



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

morst_res = MORST(G=X, obj=obj, weights= rep(1/sqrt(ngene),ngene), siglevel = 0.05)

print(paste("MORST P-Value:",morst_res))


### Measure Type-I error
expit <- function(x){
    tmp = exp(x)/(1+exp(x))
    return(exp(x)/(1+exp(x)))
}

typeI_test <- function(gi.mat, n.cnv, lavg.ln.overallmean,siglevel, no_ccret=FALSE){
    ngene       = ncol(gi.mat)
    nsubj       = nrow(gi.mat)
    idx = sample.int(nsubj, nsubj, replace=TRUE)
    X_sample = gi.mat[idx,]
    ncnv_sample = n.cnv[idx]; size_sample = lavg.ln.overallmean[idx]
    Z = cbind(ncnv_sample, size_sample)
    Z_adj = (Z - mean(Z))/sd(Z)
    beta_cnv = beta_size = log(1.5)
    case_prob = expit(-2.5 + beta_cnv*Z_adj[,1] + beta_size*Z_adj[,1])
    Y = rbinom(nsubj,1, case_prob)

    ### Run CCRET
    p_CCRET = NULL
    if(!(no_ccret)){
        GC1mat = getK.gi.mat.fun(X_sample, kernel=1, WT=rep(1/ngene, ngene))
        p_CCRET = vctest.btqt.Gmain.fun(y=Y, geno=X_sample, x.adj=Z, trait.type="binomial", SSS=GC1mat)[1]
    }
    ### Run MORST
    obj = MORST_NULL_Model(Y, Z)
    p_MORST = MORST(G=Matrix(as.matrix(X_sample),sparse=TRUE), obj=obj, weights= rep(1/sqrt(ngene),ngene), siglevel = siglevel)
    return(list(p_ccret = p_CCRET, p_morst = p_MORST))
}

typeI_df = matrix(0, nc=2, nr=3)
sigs = c(0.05, 1e-2,1e-3)
nsims = 1e4
p_ccrets = rep(NA, nsims)
for(i in 1:nrow(typeI_df)){
    cat("========================\n")
    cat("Test size:",sigs[i],"\n")
    cat("=======================\n")
    p_tmp = matrix(NA, nc=2, nr=nsims)
    sig_level = sigs[i]
    if(sig_level != 0.05) {no_ccret = TRUE} else{no_ccret=FALSE}
    for(sim in 1:nsims){
        if(sim == 1 || sim %% 100 == 0){
            cat("Simulation:",sim,"\n")
        }
        pvals = typeI_test(gi.mat, n.cnv, lavg.ln.overallmean,sig_level, no_ccret=no_ccret)
        if(no_ccret == FALSE) {
            p_tmp[sim,1] = pvals$p_ccret
            p_ccrets[sim] = pvals$p_ccret
        } 
        p_tmp[sim,2] = pvals$p_morst
    }
    p_tmp[,1] = p_ccrets
    typeI_df[i,] = c(sum(p_tmp[,1] < sig_level)/nsims, sum(p_tmp[,2] < sig_level/nsims))
}

print(typeI_df)

save(typeI_df, file="typeI_df.rds")

## Measure Power:
ngene       = ncol(gi.mat)
nsubj       = nrow(gi.mat)
idx = sample.int(nsubj, nsubj, replace=TRUE)
X_sample = gi.mat[idx,]
ncnv_sample = n.cnv[idx]; size_sample = lavg.ln.overallmean[idx]
Z = cbind(ncnv_sample, size_sample)
Z_adj = (Z - mean(Z))/sd(Z)
beta_cnv = beta_size = log(1.5)
case_prob = expit(-2.5 + beta_cnv*Z_adj[,1] + beta_size*Z_adj[,1])
Y = rbinom(nsubj,1, case_prob)

# library(foreach)
# library(parallel)
# library(doParallel)
# # Sys.setenv(TAR = "/bin/tar")
# #devtools::install_github("r-pkg-examples/rcpp-and-doparallel")
# library("Rcpp2doParallel")
# cl <- makeCluster(detectCores()-1, type="FORK")
# registerDoParallel(cl)
# res <- foreach(x=1:nsims, .combine="rbind", .packages=c("Rcpp","Rcpp")) %dopar% 
#     {print(x)
#     sourceCpp("Saddle.cpp")
#     source("Association_tests.R")
#     tmp=typeI_test(gi.mat, n.cnv, lavg.ln.overallmean,sig_level)
#     c(tmp$p_ccret, tmp$p_morst)}

# stopImplicitCluster()