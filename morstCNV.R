## MORST CNV functions and script

## Maintained by Eduardo Maury (eduardo_maury@hms.harvard.edu)
suppressMessages(library(expm))
suppressMessages(library("MASS"))
suppressMessages(library("matrixStats"))
suppressMessages(library("Matrix"))
suppressMessages(library("Rcpp"))
suppressMessages(library("RcppEigen"))
# sourceCpp("Saddle.cpp")
# source("utils.R")


old_morst_cnv <- function(Y,Z,G,kernel="dosage",weights.beta=c(1,25),weights=NULL,tau="Minimax",power.wanted=0.5,siglevel=1e-04){
    ## Main function to run morst cnv model
    n = length(Y)

    # Calculate null model
    # TODO: require fitting the null model before since we only need to fit it once
    obj <- fit_null(Y,Z)
    out_type = obj[["out_type"]]
    if(out_type=="C"){
        X.med <- obj[["X.med"]]
        Y.res <- obj[["Y.res"]]/sqrt(obj[["sigma2"]])
        n <- length(Y)
    } else if (out_type =="D"){
        X.tilde = obj[["X.tilde"]]
        Y.res = obj[["Y.res"]]
        sigma2.Y = obj[["sigma2.Y"]]
    }

    # MAF for CNVs
    MAF <- colSums(G !=2)/(nrow(G))
    p <- length(MAF)

    # weights
    if (kernel == "linear.weighted"){
        if(is.null(weights)){
            W <- dbeta(MAF, weights.beta[1],weights.beta[2])
        } else{
            if (length(weights)==p){
                W <- weights
            } else {
                stop("The length of the weights must be equal to the number of CNVRs.")
            }
        }
    } else if(kernel == "linear"){
        W <- rep(1,p)
    } else if(kernel == "dosage") {
        Kds = tcrossprod(ds.to.dsnew(X))
        E <- eigen(Kds, symmetric=TRUE)
        V <- E$vectors; U <- solve(V)
        D <- diag(E$values)
        W <- V %*% D^(1/2) %*% U
    } else{
        stop("The kernel name is not valid!")
    }

    # Allow for use sparse Matrix
    if (class(G)=="matrix"){
        if(class(W) == "matrix"){
            WG = t(G) %*%W
        } else{
            WG <- (t(G)*W)
        }
        if(out_type == "C"){
            S = as.vector(WG%*%Y.res) # score when rho=0
            Sigma0 = tcrossprod(WG) - tcrossprod(WG%*%X.med)
        } else if ( out_type=="D"){
            S = as.vector(WG%*%Y.res)
            Sigma0 = tcrossprod(t(t(WG)*sqrt(sigma2.Y)))-(WG%*%(X.tilde*sigma2.Y))%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t(WG%*%(X.tilde*sigma2.Y))
        }

    } else if(class(G)=='dgMatrix'){
        if(class(W) == "matrix"){
            GW = G%*%W
        } else{
            GW <- G%*%Diagonal(p,W)
        }
        if(out_type=="C"){
            S = as.vector(Y.res%*% GW)
            Sigma0<-(as.matrix(crossprod(GW))-crossprod(as(t(X.med)%*%GW,"matrix")))
        } else if(out_type=="D"){
            S = as.vector(Y.res %*% GW)
            tGW_X.tilde_sigma2<-as.matrix(crossprod(GW,X.tilde*sigma2.Y))
            Sigma0<-as.matrix(crossprod(GW*sqrt(sigma2.Y)))-tGW_X.tilde_sigma2%*%solve(t(X.tilde)%*%(X.tilde*sigma2.Y))%*%t(tGW_X.tilde_sigma2)
        }
    } else {
        stop("The class of genetic matrix must be a matrix or a dgCMatrix!")
    }

    # For numerical stability
    # const <- length(S)/sum(diag(Sigma0))
    # S = S/sqrt(const)
    # Sigma0 <- Sigma0*const

    ## TODO: add different rho correlation. 

    if(tau == 0){ ## SKAT
        Q <- sum(S^2)
        w.umpa <- eigen(Sigma0, symmetric=TRUE, only.values=TRUE)[["values"]]
    } else{
        eg.values <- try(eigen(Sigma0, symmetric =TRUE), TRUE)
        if(class(eg.values)=="try-error"){
            S <- S*sqrt(n)
            Sigma <- Sigma0 *n
            eg.values <- eigen(Sigma, symmetric = TRUE)
        }
        eg.vector <- eg.values[["vectors"]]
        eg.values <- eg.values[["values"]]

        eg.values[eg.values <1e-16] <- 0
        k <- max(which(eg.values>0))

        # standarize PCs
        V <- (t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%S

        ## Choose tau. If tau=0, it is SKAT.
        if (tau == "Minimax.approx"){
            tau0 = tau_ump(eg.values, power.wanted, siglevel)
        } else if (tau=="Minimax"){
            tau0 = Minmax.tau(eg.values, siglevel)[1]
        } else if(is.numeric(tau) && length(tau)==1){
            tau0 <- tau
        }
        ### weights of PCs
        w.umpa <- eg.values/(1+eg.values*tau0)
        w.umpa <- w.umpa/sum(w.umpa)
        w.umpa <- w.umpa[1:k]
        ## Q MORST
        Q <- sum(w.umpa*V^2)
    }

    pval.MORST <- Get_Q_pval(Q,w.umpa)
    
    return(pval.MORST)

}


aMORST_cnv <- function(Y,Z,G,kernel="dosage",weights.beta=c(1,25),weights=NULL,tau="Minimax",power.wanted=0.5,siglevel=1e-04){
    ## Function that perform aggregated CNV

    p_tau = morst_cnv(Y,Z,G,kernel=kernel, weights.beta=weights.beta, weights=weights, tau=tau, power.wanted=power.wanted, siglevel=siglevel)
    p_gamma = morst_cnv(Y,Z,abs(2-G),kernel=kernel, weights.beta=weights.beta, weights=weights, tau=tau, power.wanted=power.wanted, siglevel=siglevel)

    if(p_tau == 0 && p_gamma == 1){
        ## This scenario can occur when there is lack of accuracy: return tippet's method
        return(list(aMorst=min(p_tau,p_gamma), p_tau=p_tau, p_gamma=p_gamma))
    }

    if(p_tau == 1 && p_gamma == 0){
        ## This scenario can occur when there is lack of accuracy return tipett's method
        return(list(aMorst=min(p_tau,p_gamma), p_tau=p_tau, p_gamma=p_gamma))
    }

    pval.aMORST = CCT.pval(c(p_tau, p_gamma))
    return(list(aMorst=pval.aMORST, p_tau=p_tau, p_gamma=p_gamma))

}


morst_cnv <- function(Y,Z,X,null_model,kernel="dosage",weights.beta=c(1,25),weights=NULL,tau="Minimax",power.wanted=0.5,siglevel=1e-04){
    ## implementation by Eduardo Maury 2021-04-01 (april fools!)
    ## 
    n = length(Y)
    p = ncol(X)

    obj <- null_model
    out_type = obj[["out_type"]]
    if(out_type=="C"){
        X.med <- obj[["X.med"]]
        Y.res <- obj[["Y.res"]]/sqrt(obj[["sigma2"]])
        n <- length(Y)
    } else if (out_type =="D"){
        Z.tilde = obj[["X.tilde"]]
        Y.res = obj[["Y.res"]]
        sigma2.Y = obj[["sigma2.Y"]]
    }

    if(kernel =="linear.weighted"){
        MAF = colSums(X != 2)/(2*nrow(X))
        W <- diag(dbeta(MAF, 1,25))
    } else if (kernel == "linear"){
        W <- diag(p)
    } else {
        stop("The kernel name is not valid!")
    }

    ## Compute Sigma
    
    S = t(W)%*%t(X)%*%Y.res # score statistic
    H = diag(sigma2.Y) # variance function under the null. mu(1-mu) for binary response
    if(is.null(Z.tilde)){ P = diag(p)} else{
        P = H - H%*%Z.tilde%*%solve(t(Z.tilde)%*%H%*%Z.tilde)%*%t(Z.tilde)%*%H # residual generating matrix 
    }

    Sigma = t(W)%*%t(X)%*%P%*%X%*%W
    
    ## In progress: 
    ## Kernel Sigma
    # S = Y.res
    # Sigma = tcrossprod(ds.to.dsnew(X))

    ## compute the p-values for a range of rho values for a exchengeable matrix
    r.corr <- c(0,0.1,0.02,0.3,0.4,0.5)
    r.corr <- c(r.corr^2, 0.5,1)
    rho = 0.1
    E = matrix(rho, nrow=p, ncol=p)
    diag(E) = 1
    E <- chol(E)
    S = E%*%S
    Sigma = E%*%Sigma%*%E

    ### For the reason of numerical stability
    # const <- length(S)/sum(diag(Sigma))
    # S <- S*sqrt(const)
    # Sigma <- Sigma*const

    if(tau ==0){
        ## SKAT
        Q <- sum(S^2)
        w.umpa <- eigen(Sigma, symmetric = TRUE, only.values = TRUE)[["values"]]
    } else{
        ### get eigen values
        eg.values <- try(eigen(Sigma, symmetric = TRUE), TRUE)
        if (class(eg.values)=="try-error"){
            S  = S*sqrt(n)
            Sigma = Sigma*n
            eg.values = eigen(Sigma, symmetric=TRUE)
        }

        eg.vector <- eg.values[["vectors"]]
        eg.values <- eg.values[["values"]]
        eg.values[eg.values<1e-16] <- 0 
        k <- max(which(eg.values>0))

        ### The Q statistic for MORST can also be represented in a standardized PC format
        ### See section 2.7 of the Liu, Li, Lin 2020 MORST paper in JASA
        V <- (t(eg.vector[,1:k,drop=FALSE])/sqrt(eg.values[1:k]))%*%S ## V_i = t(U_i)%*%S/lambda_i^0.5

        ## perform bisection search to find minimax tau
        tau0 = tau_ump(eg.values, power.wanted, siglevel=siglevel)

        ### Caclulate the weights lambda/(1+lambda*tau0) for the chi square stat
        w.umpa <- eg.values/(1+ eg.values*tau0)
        w.umpa <- w.umpa/sum(w.umpa)
        w.umpa <- w.umpa[1:k]
        Q <- sum(w.umpa*V^2)
    }

    pval.MORST <- Get_Q_pval(Q,w.umpa)



}
