

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

Get.alpha0<-function(logOR,prevalance){
    ###Get the upper and lower values for alpha0, then bisectional search
    alpha0<-log(prevalance/(1-prevalance))-mean(logOR)
    preva<-mean(1-1/(1+exp(logOR+alpha0)))
    if (preva>=prevalance){
        alpha0.upper<-alpha0
        a<-1
        while (preva>=prevalance){
            alpha0<-log(prevalance/(1-prevalance))-mean(logOR)-a*abs(mean(logOR))
            preva<-mean(1-1/(1+exp(logOR+alpha0)))
            a<-a+1
        }
        alpha0.lower<-alpha0
    }else{
        alpha0.lower<-alpha0
        a<-1
        while (preva<prevalance){
            alpha0<-log(prevalance/(1-prevalance))-mean(logOR)+a*abs(mean(logOR))
            preva<-mean(1-1/(1+exp(logOR+alpha0)))
            a<-a+1
        }
        alpha0.upper<-alpha0
    }
    #### bisectional search
    acc<-prevalance*0.01
    while (abs(preva-prevalance)>acc){
        alpha0<-(alpha0.lower+alpha0.upper)/2
        preva<-mean(1-1/(1+exp(logOR+alpha0)))
        if (preva>=prevalance){
            alpha0.upper<-alpha0
        }else{
            alpha0.lower<-alpha0
        }
    }
    return(alpha0)
}




MORST_NULL_Model<-function(Y,X=NULL){
    n<-length(Y)
    #### check the type of Y
    if ((sum(Y==0)+sum(Y==1))==n){
        out_type<-"D"
    }else{
        out_type<-"C"
    }
    #### Add intercept
    X.tilde<-cbind(rep(1,length(Y)),X)
    #colnames(X.tilde)[1]<-"Intercept"
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

MORST<-function(G,obj,kernel="linear.weighted",weights.beta=c(1,25),r.corr=0,weights=NULL,tau="Minimax.approx",power.wanted=0.5,siglevel=1e-04){
    ### check obj
    if (names(obj)[1]!="out_type"){
        stop("obj is not calculated from MORST_NULL_MODEL!")
    }else{
        out_type<-obj[["out_type"]]
        if (out_type=="C"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.med","Y.res","sigma2"))){
                stop("obj is not calculated from MORST_NULL_MODEL!")
            }else{
                X.med<-obj[["X.med"]]
                Y.res<-obj[["Y.res"]]/sqrt(obj[["sigma2"]])  ## rescaled residules
                n<-length(Y)
            }
        }else if (out_type=="D"){
            if (!all.equal(names(obj)[2:length(obj)],c("X.tilde","Y.res","sigma2.Y"))){
                stop("obj is not calculated from MORST_NULL_MODEL!")
            }else{
                X.tilde<-obj[["X.tilde"]]
                Y.res<-obj[["Y.res"]]
                sigma2.Y<-obj[["sigma2.Y"]]
                n<-length(Y.res)
            }
        }
    }
    ### MAF
    MAF<-colSums(G)/(2*dim(G)[1])
    p<-length(MAF)
    #### weights
    if (kernel=="linear.weighted"){
        if (is.null(weights)){
            W<-dbeta(MAF,weights.beta[1],weights.beta[2])
        }else{
            if (length(weights)==p){
                W<-weights
            }else{
                stop("The length of weights must equal to the number of variants!")
            }
        }
        
    }else if (kernel=="linear"){
        W<-rep(1,p)
    }else{
        stop("The kernel name is not valid!")
    }
    ###### if G is sparse or not
    if (class(G)=="matrix"){
        #### WG
        WG<-(t(G)*W)   # diag(W)%*% t(G)
        #### The z statsitics and the covariance Sigma
        if (out_type=="C"){
            ### Z statistics when rho=0
            Z.stat0<-as.vector(WG%*%Y.res)
            ### Sigma when rho=0
            Sigma0<-(tcrossprod(WG)-tcrossprod(WG%*%X.med))
        }else if (out_type=="D"){
            ### Z statistics when rho=0
            Z.stat0<-as.vector(WG%*%Y.res)
            ### Sigma when rho=0
            Sigma0<-tcrossprod(t(t(WG)*sqrt(sigma2.Y)))-(WG%*%(X.tilde*sigma2.Y))%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t(WG%*%(X.tilde*sigma2.Y))
        }
    }else if (class(G)=="dgCMatrix"){
        #### GW
        W.diag<-Diagonal(p,W)
        GW<-(G%*%W.diag)
        #### The z statsitics and the covariance Sigma
        if (out_type=="C"){
            ### Z statistics when rho=0
            Z.stat0<-as.vector(Y.res%*%GW)
            ### Sigma when rho=0
            Sigma0<-(as.matrix(crossprod(GW))-crossprod(as(t(X.med)%*%GW,"matrix")))
        }else if (out_type=="D"){
            ### Z statistics when rho=0
            Z.stat0<-as.vector(Y.res%*%GW)
            ### Sigma when rho=0
            tGW_X.tilde_sigma2<-as.matrix(crossprod(GW,X.tilde*sigma2.Y))
            Sigma0<-as.matrix(crossprod(GW*sqrt(sigma2.Y)))-tGW_X.tilde_sigma2%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t(tGW_X.tilde_sigma2)
        }
    }else{
        stop("The class of G must be matrix or dgCMatrix!")
    }
    
    ### For the reason of numerical stability
    const<-length(Z.stat0)/sum(diag(Sigma0))
    Z.stat0<-Z.stat0*sqrt(const)
    Sigma0<-Sigma0*const
    
    
    ### compute the p-values for a range of rho values 
    if (is.null(r.corr)){
        r.corr<-c(0,0.1,0.2,0.3,0.4,0.5)
        r.corr<-c(r.corr^2,0.5,1)
    }
    
    pvals<-c()
    for (rho in r.corr){
        if (rho > 1 || rho < 0){
            stop("Values in r.corr should be between 0 and 1!")
        }else if (rho == 1){
            V<-sum(Z.stat0)/sqrt(sum(Sigma0))
            w.umpa<-1 ## degree of freedom
            Q<-V^2   ## Q test statistic
        }else{
            if (rho==0){
                Z.stat<-Z.stat0 ### Z statistics
                Sigma<-Sigma0   ### Sigma
            }else{
                E<-matrix(rho,nrow = p,ncol = p)
                diag(E)<-1
                E<-chol(E)
                ### Z statistics
                Z.stat<-E%*%Z.stat0
                ### Sigma
                Sigma<-E%*%Sigma0%*%t(E)
            }
            
            if (tau==0){ ## SKAT
                Q<-sum(Z.stat0^2)
                w.umpa<-eigen(Sigma0,symmetric = TRUE,only.values = TRUE)[["values"]]
            }else{
                ### eigenvalue decomposition 
                eg.values<-try(eigen(Sigma,symmetric = TRUE),TRUE)
                if (class(eg.values)=="try-error"){
                    Z.stat<-Z.stat*sqrt(n)
                    Sigma<-Sigma*n
                    eg.values<-eigen(Sigma,symmetric = TRUE)
                }
                eg.vector<-eg.values[["vectors"]]
                eg.values<-eg.values[["values"]]
                
                eg.values[eg.values<1e-16]<-0
                
                k<-max(which(eg.values>0))
                
                ##### standarized PCs 
                V<-(t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%Z.stat
                ##### choose tau. If tau=0, it is SKAT.
                if (tau=="Minimax.approx"){
                    tau0<-tau_ump(eg.values,power.wanted,siglevel)
                }else if (tau=="Minimax"){
                    tau0<-Minmax.tau(eg.values,siglevel)[1]
                }else if (is.numeric(tau) && length(tau)==1){
                    tau0<-tau
                }
                #### weights for PCs
                w.umpa<-eg.values/(1+eg.values*tau0)
                w.umpa<-w.umpa/sum(w.umpa)
                w.umpa<-w.umpa[1:k]
                #### Q
                Q<-sum(w.umpa*V^2)
            }
            
        }
        pval.MORST<-Get_Q_pval(Q,w.umpa)
        pvals<-c(pvals,pval.MORST)
    }
    ##### use the Cauchy method to combine p-values
    if (length(pvals)>1){
        pval.MORST.O<-CCT.pval(pvals)
    }else{
        pval.MORST.O<-pvals
    }
    
    return(pval.MORST.O)
}




CheckWeights<-function(weights,p){
    if (is.vector(weights)){
        if (length(weights)!=p){
            stop("The length of weights must be equal to number of variants!")
        }else{
            weights<-as.matrix(weights,ncol=1)
        }
    }
    
    if (!is.matrix(weights)){
        stop("weights must be a matrix or vector!")
    }else if (nrow(weights)!=p){
        stop("The number of rows in weights must be equal to number of variants!")
    }
    
    if (mode(weights)!="numeric"){
        stop("weights must be numeric!")
    }else if (sum(weights<0)>0){
        stop("All the weights should be non-nagetive!")
    }else if (sum(colSums(weights)==0)>0){
        stop("At least one weight must be positive!")
    }
    
    ### standarize weights
    weights.standarized<-weights
    for (k in 1:ncol(weights)){
        weights.standarized[,k]<-weights.standarized[,k]/sum(weights[,k])*p
    }
    
    return(weights.standarized)
}

