
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
    typeI_df[i,] = c(sum(p_tmp[,1] < sig_level)/nsims, sum(p_tmp[,2] < sig_level)/nsims))
}

print(typeI_df)

save(typeI_df, file="typeI_df.rds")

#####################################################################
#######################################################################
## Measure Power:

get_dupdelonly_genes <- function(df){
    duponly_genes = c()
    delonly_genes = c()
    for(i in ncol(df)){
        if(all(df[,i] == 3 | df[,i] == 0)){
            duponly_genes = c(duponly_genes,names(df[,i,drop=FALSE]))
        } else if(all(df[,i] == 2 | df[,i] == 0)){
            delonly_genes = c(delonly_genes, names(df[,i,drop=FALSE]))
        }
    }
    return(list(duponly_genes=duponly_genes, delonly_genes=delonly_genes))
}

only_genes = get_dupdelonly_genes(gi.mat)


power_test <- function(effect_size,gi.mat, n.cnv, lavg.ln.overallmean,siglevel, no_ccret=FALSE){
    ngene       = ncol(gi.mat)
    nsubj       = nrow(gi.mat)
    num_causative = ceiling(ngene*0.05)
    idx = sample.int(nsubj, nsubj, replace=TRUE)
    X_sample = gi.mat[idx,]
    ncnv_sample = n.cnv[idx]; size_sample = lavg.ln.overallmean[idx]
    Z = cbind(ncnv_sample, size_sample)
    Z_adj = (Z - mean(Z))/sd(Z)
    alpha_cnv = alpha_size = log(1.5)
    beta_gi = matrix(c(rep(effect_size, num_causative), rep(0, ngene - num_causative)), nc=1)
    case_prob = expit(-2.5 + alpha_cnv*Z_adj[,1] + alpha_size*Z_adj[,1] + as.matrix(X_sample)%*%beta_gi)
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
    #return(list(p_ccret = p_CCRET[1], p_morst = p_MORST))
    return(c(as.numeric(p_CCRET),p_MORST))
}

power_range = seq(0.5,1.5, length.out=5)
power_df = matrix(0, nc=3, nr=length(power_range))
sigs = c(0.05, 1e-4)
nsims = 300
p_ccrets = rep(NA, nsims)
gi.mat[gi.mat > 0] = 1
library(purrr)
sig_level = 1e-4
for(i in 1:nrow(power_df)){
    effect_size = power_range[i]
    cat("========================\n")
    cat("Effect size:", effect_size,"\n")
    cat("=======================\n")
    p_tmp = matrix(NA, nc=2, nr=nsims)
    # for(sim in 1:nsims){
    #     if(sim == 1 || sim %% 100 == 0){
    #         cat("Simulation:",sim,"\n")
    #     }
    #     pvals = power_test(effect_size,gi.mat, n.cnv, lavg.ln.overallmean,sig_level)
    #     p_tmp[sim,1] = pvals$p_ccret
    #     p_tmp[sim,2] = pvals$p_morst
    # }
    p_tmp = nsims %>% rerun(power_test(effect_size,gi.mat, n.cnv, lavg.ln.overallmean,sig_level)) %>% reduce(rbind)
    power_df[i,] = c(sum(p_tmp[,1] < sig_level)/nsims, sum(p_tmp[,2] < sig_level)/nsims, sig_level)        
}


# ## attempt with mcapply
# library(parallel)
# cl <- makeCluster(detectCores()-1, type="FORK")


# system.time(10 %>% rerun(power_test(effect_size,gi.mat, n.cnv, lavg.ln.overallmean,sig_level)) %>% reduce(rbind))
# system.time( for (i in 1:10){power_test(effect_size,gi.mat, n.cnv, lavg.ln.overallmean,sig_level)})

# ## alpha < 1e-4, singal 0.4
# > power_df
#      [,1]  [,2]  [,3]
# [1,] 0.17 0.245 1e-04
# [2,] 0.52 0.640 1e-04
# [3,] 0.87 0.940 1e-04
# [4,] 0.96 1.000 1e-04
# [5,] 1.00 1.000 1e-04



## plotting power curve
power_df %>% pivot_longer(cols=c(1,2), names_to="model", values_to="Power") %>% 
    ggplot(aes(x=beta, y=Power)) + 
    geom_point(aes(color=model, shape=model), size=2) + 
    geom_line(aes(color=model), lwd=1) + theme_bw(base_size=12) + 
    scale_color_manual(values=c("CCRET" = "orange","CNV.MORST" = "dodgerblue3")) + 
    labs(title = "test-size: 1e-4")

###################################################################################
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