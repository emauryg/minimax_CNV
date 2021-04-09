## Code for power simulations


## Maintained by Eduardo Maury (eduardo_maury@hms.harvard.edu)
suppressMessages(library("MASS"))
suppressMessages(library("matrixStats"))
suppressMessages(library("Matrix"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppEigen"))

# sourceCpp("Saddle.cpp")
# source("utils.R")
# source("morstCNV.R")


## input:
##      calls: a data.frame with columuns (chr, start, end, sample_id, score)

library(parallel)

power_sim_parallel <- function(i,Xmat,idx_names, num_loci=100, dup_effect=0.5,del_effect=-0.5, IBS=FALSE){
    ## Generate random samples
    draw_names <- sample(idx_names, size=2000, replace=FALSE)
    X_sim = Xmat[draw_names,]
    #reduce to areas that have at least one event
    X_sim = X_sim[,1:num_loci]
    actual.freq = colMeans(X_sim != 2)
    X_sim = X_sim[,actual.freq  >0]
    n = nrow(X_sim)
    p = ncol(X_sim)
    Y = rep(0,n)
    signal_frac = as.integer(floor(0.2*p))
    ## Simulate X,Y,Z
    Z<-cbind(rnorm(n),(rbinom(n,size = 1,prob = 0.5)-0.5))
    loc.signal <- sample.int(p, signal_frac)
    beta_dup = dup_effect*(rbinom(length(loc.signal),size=1, prob=0.5))
    beta_del = del_effect*(rbinom(length(loc.signal),size=1, prob=0.5))
    case_prob = expit(Z[,1]*0.5 + Z[,2]*0.5 + (X_sim[,loc.signal, drop=FALSE] == 3) %*% beta_dup + 
                      (X_sim[,loc.signal,drop=FALSE]==1) %*% beta_del)
    Y = rbinom(length(case_prob), size=1, prob=case_prob)
    Kmat.ds = tcrossprod(ds.to.dsnew(X_sim))
    null_model = fit_null(Y,Z=Z)
    res1 = aMORST_cnv(Y,Z=Z, X_sim, null_model,IBS=IBS, kernel="linear", tau="Minimax.approx", power.wanted=0.5, siglevel=0.05)
    res2 = vctest.btqt.Gmain.fun(y=Y, geno=X_sim, x.adj =Z, trait.type="binomial", SSS=Kmat.ds)
    #res2 = 0
    pvals = c(res2,res1$aMorst, res1$p_tau, res1$p_gamma)
    #message(paste(pvals[i,],sep = "\t"))
    return(pvals)
}

do_power <- function(X,nsims, num_loci=100, dup_effectRange, del_effectRange){
    if(length(dup_effectRange) != length(del_effectRange)){
        stop("The ranges of gains and losses must be of the same size!")
    }
    n_iter = length(dup_effectRange)
    power_table = matrix(0, nr=n_iter, nc=4)
    cont_names = rownames(X)
    for(i in 1:n_iter){
        message("Iteration:", i)
        tmp = mclapply(1:nsims,power_sim_parallel,Xmat=X, idx_names=cont_names,num_loci = 100,
                                   dup_effect = 6,del_effect = 6,IBS=FALSE,mc.cores = 8)
        tmp2 = do.call(rbind, tmp)
        power_table[i,] = colMeans(tmp2 < 0.05)/nsims
    }

    return(power_table)

}

## Power simulation with PGC CNV control samples

suppressMessages(library(tidyverse))
suppressMessages(library(CNVRanger))
## Download all the packages needed to run MORST and CCRET
library("MASS")
library("CompQuadForm")
library("matrixStats")
library("Matrix")
library("Rcpp")
library("RcppEigen")


## Load the PGC dataset
pgc_cnvs <- read_tsv("PGC_SCZ_callset_cnv_liftOver.bed",col_names=c("chrom","start","end","sample_id","diagnosis","state"))

## Change encoding from G and L to 3 and 1 for dosage. 

pgc_cnvs$score = 0
pgc_cnvs$score[pgc_cnvs$state=="G"] = 3
pgc_cnvs$score[pgc_cnvs$state=="L"] = 1
pgc_cnvs$state[pgc_cnvs$state=="G"] = "gain"
pgc_cnvs$state[pgc_cnvs$state=="L"] = "loss"

pgc_cnvs = pgc_cnvs %>% mutate(batch = str_split(sample_id,pattern = "[*]",simplify = TRUE)[,1]) %>%
    mutate(batch = str_split(batch, "_",simplify=TRUE)[,5])

calls = pgc_cnvs %>% dplyr::select(c(chrom,start, end, sample_id,score))
colnames(calls)[1:3] = c("chr","start","end")

cont_names = pgc_cnvs$sample_id[pgc_cnvs$diagnosis=="CONT"]
calls = calls[match(c(cont_names), calls$sample_id),]
grl = makeGRangesListFromDataFrame(calls, split.field="sample_id", keep.extra.columns=TRUE)
ra =  RaggedExperiment(grl)
X = disjoinAssay(ra, simplifyDisjoin=mean)
X[is.na(X)] = 2
X = t(X)

power_res = do_power(X, 2000, num_loci=100,
                    dup_effectRange = c(0.1,0.5,0.75,1),
                    del_effectRange = c(0.1,0.5,0.75,1))



