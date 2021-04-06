## Utility functions for morst CNV

fit_null <- function(Y,Z=NULL){
    ## Function to estimate \hat{mu_0}
    ## input:
    ##      Y: outcome vector (can onlyhandle continuous and binary for now)
    ##      Z: matrix of covariates that we want to control for, dim n X p
    
    n = length(Y)
    # Check the type of Y
    if ((sum(Y==0)+sum(Y==1))== n){
        out_type = "D"
    } else{
        out_type = "C"
    }

    # Add intercept of the model does not have one already
    if (!is.null(Z) && all(Z[,1] == 1)){
        X.tilde = Z
    } else{
        X.tilde = cbind(rep(1,n), Z)
    }

    if (out_type=="C"){
        #### estimate of sigma square 
        X.med<-X.tilde%*%solve(chol(t(X.tilde)%*%X.tilde))   ## X.med%*%t(X.med) is the projection matrix of X.tilde
        Y.res<-as.vector(Y-(Y%*%X.med)%*%t(X.med))
        sigma2<-sum(Y.res^2)/(n-ncol(X.med))
        #### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["X.med"]]<-X.med
        res[["Y.res"]]<-Y.res
        res[["sigma2"]]<-sigma2
    }else if (out_type=="D"){
        #### fit null model
        g<-glm(Y~0+X.tilde,family = "binomial")
        prob.est<-g[["fitted.values"]]
        #### unstandarized residuals
        Y.res<-(Y-prob.est)
        ### Sigma when rho=0
        sigma2.Y<-prob.est*(1-prob.est)  ### variance of each Y_i
        ### output
        res<-list()
        res[["out_type"]]<-out_type
        res[["X.tilde"]]<-X.tilde
        res[["Y.res"]]<-Y.res
        res[["sigma2.Y"]]<-sigma2.Y
    }
    return(res)
}




Get_Q_pval<-function(Q,w){
    ############## when length(w)==1
    if (length(w)==1){
        Q<-Q/w
        w<-1
        pval<-1-pchisq(Q,df=w)
        if (pval<=0){
            pval<-liu(Q,w)
        }
        return(pval)
    }
    
    ############### when length(w)>1
    lim.davies.up<-1e+09
    #### when pval>1e-05, try davies method
    pval<-suppressWarnings(davies(Q,w,acc=1e-06))
    if (pval[["ifault"]]!=0){
        pval<-davies(Q,w,lim=lim.davies.up,acc=1e-06)
        if (pval[["ifault"]]!=0){
            stop("Davies method has an error!")
        }else{
            pval<-pval[["Qq"]]
        }
    }else{
        pval<-pval[["Qq"]]
    }
    # if good, then return
    if (pval>1e-05){
        return(pval)
    }
    
    ##### when 1e-10>pval>1e-05, try davies method
    pval<-suppressWarnings(davies(Q,w,lim=1e+06,acc=1e-11))
    if (pval[["ifault"]]!=0){
        pval<-davies(Q,w,lim=lim.davies.up,acc=1e-11)
        if (pval[["ifault"]]!=0){
            stop("Davies method has an error!")
        }else{
            pval<-pval[["Qq"]]
        }
    }else{
        pval<-pval[["Qq"]]
    }
    # if good, then return
    if (pval>1e-10){
        return(pval)
    }
    ##### when 1e-15>pval>1e-10, try saddle point method
    pval<-Saddle(Q,w)
    # if good, then return
    if (pval>1e-15){
        return(pval)
    }
    ##### when pval<1e-15 use liu's method
    pval<-liu(Q,w)
    return(pval)
}


CCT.pval<-function(Pvals,Weights=NULL){
    #### check if there is NA
    if (sum(is.na(Pvals))>0){
        stop("Cannot have NAs in the p-values!")
    }
    #### check if Pvals are between 0 and 1
    if ((sum(Pvals<0)+sum(Pvals>1))>0){
        stop("P-values must be between 0 and 1!")
    }
    #### check if there are pvals that are either exactly 0 or 1.
    is.zero<-(sum(Pvals==0)>=1)
    is.one<-(sum(Pvals==1)>=1)
    if (is.zero && is.one){
        stop("Cannot have both 0 and 1 p-values!")
    }
    if (is.zero){
        return(0)
    }
    if (is.one){
        warning("There are p-values that are exactly 1!")
        return(1)
    }
    
    #### Default: equal weights. If not, check the validity of the user supplied weights and standadize them.
    if (is.null(Weights)){
        Weights<-rep(1/length(Pvals),length(Pvals))
    }else if (length(Weights)!=length(Pvals)){
        stop("The length of weights should be the same as that of the p-values")
    }else if (sum(Weights<0)>0){
        stop("All the weights must be positive!")
    }else{
        Weights<-Weights/sum(Weights)
    }
    
    
    #### check if there are very small non-zero p values
    is.small<-(Pvals<1e-16)
    if (sum(is.small)==0){
        cct.stat<-sum(Weights*tan((0.5-Pvals)*pi))
    }else{
        cct.stat<-sum((Weights[is.small]/Pvals[is.small])/pi)
        cct.stat<-cct.stat+sum(Weights[!is.small]*tan((0.5-Pvals[!is.small])*pi))
    }
    #### check if the test statistic is very large.
    if (cct.stat>1e+15){
        pval<-(1/cct.stat)/pi
    }else{
        pval<-1-pcauchy(cct.stat)
    }
    return(pval)
}


expit <- function(x){
    return(exp(x)/(1+exp(x)))
}