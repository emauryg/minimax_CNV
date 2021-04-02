########################################################################################
## (2)  FUNCTIONS
########################################################################################
library(MASS)
library(CompQuadForm)



##============================================
##   get kernel matrix
##============================================
ds.to.dsnew<-function(ds)
{
  ncols = ncol(ds)
  nrows = nrow(ds)
  dsnew = matrix(0, ncol=(3 * ncols), nrow=nrows)
  for(j in 1:ncols) 
  {
    jj = 3*(j-1)+1
    dsnew[ds[,j] <2, jj  ] = 1
    dsnew[ds[,j]==2, jj+1] = 1
    dsnew[ds[,j] >2, jj+2] = 1
  }
  return(dsnew)
}


getK.gi.mat.fun=function(mat = gi.all.mat, kernel = 1, WT=rep(1/nloci, nloci))
  {
        ## mat = design mat, nsubj by nloci
        ## kernel =1 : linear kernel
        ## kenrel =2 : qudratic kernel
    Kroot = t(t(mat) * sqrt(WT))
    Km    = tcrossprod(Kroot)      
    ## print(range(Km - mat %*% diag(WT) %*% t(mat)))
    
    if(kernel==2){      Km = (1+ Km)^2 -1   }
    return("Km.gi"=Km)
  }




##============================================
##    VC test
##============================================
vctest.btqt.Gmain.fun<-function(y=yy, geno=Geno, x.adj=xx, trait.type="gaussian", SSS=SS3)
     {
      trait.int <- charmatch(trait.type, c("gaussian", "binomial"))
      if (is.na(trait.int))        stop("Invalid trait type")
      if (trait.int == 0)          stop("Ambiguous trait type")
      if (length(y) != nrow(geno)) stop("Dims of y and geno are not compatible")
      
      n.subj <- length(y)
      n.ij   <- n.subj *(n.subj-1)/2
      n.loci <- ncol(geno)

      adjusted <- TRUE ;        if (all(is.na(x.adj)))       adjusted <- FALSE
      if (adjusted)
        {
          ## ## standardized X
          ## x.adj = (x.adj-mean(x.adj))/sd(x.adj)
          x.adj <- as.matrix(x.adj)
          if (nrow(x.adj) != length(y)) stop("Dims of y and x.adj are not compatible")
        }
       miss <- is.na(y)
       if (adjusted)       miss <- miss | apply(is.na(x.adj), 1, any)
       y    <- as.numeric(y[!miss])
       geno <- geno[!miss, ]
       if (adjusted)       x.adj <- x.adj[!miss, , drop = FALSE]
      ##-----------------------------
      ##-----------------------------
      ##-----------------------------
      if (!adjusted){ newx.adj = as.matrix(rep(1, n.subj))      }
      if (adjusted){  newx.adj = as.matrix( cbind(rep(1, n.subj), x.adj) )}
      ##########################################
      ## marginal G test
      ##########################################
      ## source("/Users/jytzeng/Research/HSregression/Software/hsreg.fun.r")
      ## source("/Users/jytzeng/Research/HSregression/Software/hsreg.conversion.r")
      source("CCRET/hsreg.fun.r")
      source("CCRET/hsreg.conversion.r")
      ## geno = hsreg.prepare.geno.data(snp.to.impute(geno), "impute")
      ## tmp = hsreg.pval.fun(y, geno, x=x.adj, trait.type, approx.method=3)
      ## tmp = hsreg.pval.jy.fun(y, geno, x=x.adj, trait.type, approx.method=3)
      ##---------------------------
      pmat.data     = hsreg.pmat.fun(y=y, x=x.adj, trait.type=trait.type)
      P.mat       = pmat.data$pmat
      ## P.mat.class = pmat.data$pmat.class
      newy          = pmat.data$newy
      Vy            = pmat.data$vy ## JY change $v to $vy
      W  = pmat.data$w ## JY change $v to $vy
      ## --- W is a scaler if no covariates (regardless BT or QT) or if QT (regardless w/ or w/o covariate);
      ## --- W is a vector if BT with covaraites, in that case CC0's calculation need to be careful
      T.G = 1/2 * newy %*% SSS %*% newy;
      PS          = crossprod(P.mat, SSS)  ## for PS= P.mat %*% SSS ## <--CSCS originally "PS= ps.calc.2(newgeno, P.mat)
      PSP = (PS %*% P.mat)
      ## if(trait.int == 2) ## binomial
      ##   {
      ##     CC0 = 1/2 / Vy * PSP
      ##   }else if(trait.int == 1){## gaussian
      ##     CC0 = 1/2 * Vy * PSP
      ##   }
      ## ##--- the above is the same as below: -----
      if(length(W)==1){ ## if no covariates or QT (i.e., W is a scaler)
        CC0 = 1/2 / W * PSP
        eg = eigen(CC0, symmetric=T, only.values=T);
        ## eg = eigen(CC0, symmetric=T);
      }else if (length(W)>1){## if BT with covariate (i.e., W is a vector)
        CC0 = 1/2 * diag(1/sqrt(W))%*% PSP %*% diag(1/sqrt(W))
        ## eg = eigen(CC0, symmetric=T, only.values=T);
        eg = eigen(CC0, symmetric=T);
      }
      evalue = eg$values; 
          le = evalue[1]
          if (le == 0) le = 1
          ev   = evalue[(evalue / le) > 1e-10]
      c1 = sum(ev);          c2 = sum(ev^2);          c3 = sum(ev^3);          dof = hprime = c2^3/c3^2
      
      ## pval.3M = 1-pchisq( (T.G-c1)*sqrt(hprime/c2)+hprime, dof)

      #pval.T.G.liu    = liu(   T.G, ev, rep(1, length(ev)), rep(0, length(ev)))
      #pval.T.G.davies = davies(T.G, ev, rep(1, length(ev)), rep(0, length(ev)))$Qq
      pval.T.G.davies = Get_Q_pval(T.G,ev) ## using the Lin's lab optimized Davies function. 
      ##-----------------
      return(c(
               ##          "T.G" = T.G,   
               ##"pmat.data"= pmat.data,
               "pval.T.G.davies" = pval.T.G.davies
               #"pval.T.G.liu" = pval.T.G.liu
               ))
      ## return("pval.T.G.davies" = pval.T.G.davies)
    }



