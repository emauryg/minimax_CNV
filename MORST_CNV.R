
## R code to run MORST_CNV, based on  great part from https://www.tandfonline.com/doi/full/10.1080/01621459.2020.1831926

## Written by Eduardo Maury (eduardo_maury AT hms DOT harvard DOT edu)

setwd("C:/Users/emaur/Dropbox (MIT)/PhD Classes/STAT364 Scalable Statistical Inference for Big Data and Applications/presentation_Lin/minimax_CNV")

library(tidyverse)

pgc_cnvs <- read_tsv("../PGC_SCZ_callset_cnv_liftOver.bed",col_names=c("chrom","start","end","sample_id","diagnosis","state"))

## Change encoding from G and L to 3 and 1 for dosage. 

pgc_cnvs$score = 0
pgc_cnvs$score[pgc_cnvs$state=="G"] = 3
pgc_cnvs$score[pgc_cnvs$state=="L"] = 1
pgc_cnvs$state[pgc_cnvs$state=="G"] = "gain"
pgc_cnvs$state[pgc_cnvs$state=="L"] = "loss"
## Getting regions. 

library(CNVRanger)


calls = pgc_cnvs %>% dplyr::select(c(chrom,start, end, sample_id,score))
colnames(calls)[1:3] = c("chr","start","end")


grl = makeGRangesListFromDataFrame(calls, split.field="sample_id", keep.extra.columns=TRUE)

#grl = sort(grl)


gr = unlist(grl)

cnvrs = disjoin(gr)

testr = cnvrs[seqnames(cnvrs)=="chr2"]
cnvr = as.data.frame(testr)

X  = matrix(2, nc=nrow(cnvr),nr=length(grl))
n = length(grl)
p = nrow(cnvr)
for (j in 1:p){
  for (i in 1:n){
    df = as.data.frame(grl[[i]])
    if(cnvr$seqnames[j] %in% df$seqnames){
      tmp = df[df$seqnames == cnvr$seqnames[j],]
      for(l in 1:nrow(tmp)){
        if(tmp$start[l] <= cnvr$start[j] & tmp$end[l] >= cnvr$end[j]){
          X[i,j] = tmp$type[l]
        }
      }
    }
  }
  
}

ra =  RaggedExperiment(grl[seqnames(grl)=="chr2"])
X = disjoinAssay(ra, simplifyDisjoin=mean)
X[is.na(X)] = 2
X = t(X)
dim(X)
#reduce to areas that have at least one event
xmeans = colMeans(X)
X = X[,xmeans != 2]
dim(X)

Y = pgc_cnvs$diagnosis[match(rownames(X),pgc_cnvs$sample_id)]
Y = ifelse(Y=="SCZ",1,0)

res = summary(glm.fit(X,Y, family=binomial("logit")))

