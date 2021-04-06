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

power_sim <- function(X,idx_names,nsims, num_loci=100, dup_effect=0.5,del_effect=-0.5){

    power_res = c(0,0,0,0)
    ## iterate through nsims iterations and compute power
    pvals = matrix(0,nr=nsims, ncol=length(power_res))
    for(i in 1:nsims){
        ## Generate random samples
        draw_names <- sample(idx_names, size=2000, replace=FALSE)
        X_sim = X[draw_names,]
        #reduce to areas that have at least one event
        X_sim = X_sim[,1:num_loci]
        actual.freq = colMeans(X_sim != 2)
        X_sim = X_sim[,actual.freq  >0]
        n = nrow(X_sim)
        p = ncol(X_sim)
        Y = rep(0,n)
        signal_frac = as.integer(floor(0.05*p))
        beta_dup = rep(dup_effect,signal_frac)
        beta_del = rep(del_effect,signal_frac)
        for(j in 1:n){
            case_prob = expit(-2 + beta_dup*sum(X_sim[j,]==3) + beta_del*sum(X_sim[j,]==1))
            Y[j] = rbinom(1,1,case_prob)
        }
        Kmat.ds = tcrossprod(ds.to.dsnew(X_sim))
        null_model = fit_null(Y,Z=NULL)
        res1 = aMORST_cnv(Y,Z=NULL, X_sim, null_model, kernel="linear.weighted", tau="Minimax.approx", power.wanted=0.8, siglevel=1e-4)
        res2 = vctest.btqt.Gmain.fun(y=Y, geno=X_sim, x.adj = rep(1,n), trait.type="binomial", SSS=Kmat.ds)
        pvals[i,] = c(res2,res1$aMorst, res1$p_tau, res1$p_gamma)
    }

    power_res = colMeans(pvals < 1e-4)/nsims
    return(power_res)

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
        power_table[i,] = power_sim(X,idx_names=cont_names, nsims, num_loci=num_loci, 
                                    dup_effect=dup_effectRange[i], del_effect=del_effectRange[i])
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



