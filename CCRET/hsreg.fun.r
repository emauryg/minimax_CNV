library(MASS)

if (!exists("hsreg.verbosity"))
  hsreg.verbosity = 0

hsreg.timings = 0
hsreg.ps.calc = 0

hsreg.banner = function()
{
  cat("*************************************\n");
  cat("HSreg R version 1.13\n");
  cat("*************************************\n");
}

##------------------------------------------------------------------------------------------------------
## Check that a valid trait type was specified.
## Returns 1 for gaussian, 2 for binomial.
##------------------------------------------------------------------------------------------------------
hsreg.check.trait <- function(trait.type)
{
  trait.int <- charmatch(trait.type, c("gaussian", "binomial"))
  if (is.na(trait.int)) stop(gettextf("invalid trait type: %s", trait.type))
  if (trait.int == 0) stop(gettextf("ambiguous trait type: %s", trait.type))
  return(trait.int)
}

##------------------------------------------------------------------------------------------------------
## Prepare the trait data by removing entries associated with missing data in y or x.
## trait.type must be "binomial" or "gaussian".
## Returns a list with named entries:
##   y         = trait (y) data with missing entries removed.
##   x         = covariate (x) data with missing entries removed.
##   missing   = list of indices of rows that have been removed from the original data.
##------------------------------------------------------------------------------------------------------
hsreg.prepare.trait.data <- function(y, x=NA, trait.type)
{
  # Check for a valid trait type. If this is "binomial" we also 
  # check that the data is represented as 0=control,1=case,NA=missing. 
  # If we find the PLINK format (0=missing,1=control,2=case) we convert.
  trait.int = hsreg.check.trait(trait.type)
  if (trait.int == 2) # Binomial
  {
    if (all(y %in% c(NA,0,1))) {

      # Data coded as NA,0,1. While this could happen in data intended to 
      # be in PLINK format (i.e. only 0 and 1 present), it would mean 
      # that there are no cases, so should not occur.
      if (hsreg.verbosity > 0)
        cat("Trait format is: 0=control, 1=case NA=missing\n")

    } else if ( all(y %in% c(NA,0,1,2)) ) {

      # PLINK style trait values. Convert to NA,0,1.
      y[y==0] = NA
      y = y - 1

      if (hsreg.verbosity > 0)
        cat("Trait format is: 1=control, 2=case, 0=missing\n")

    } else {
      uy = sort(unique(y), na.last=TRUE)
      stop(gettextf("Unexpected values found in binomial trait data: %s", paste(uy, collapse=" ")))
    }
  }

  # Check for covariates, and whether the number of covariates matches the number of individuals.
  adjusted = !all(is.na(x))
  if (adjusted)
  {
    x = as.matrix(x)
    if (nrow(x) != length(y)) 
      stop(gettextf("Dimensions of trait data (%d) and covariate data (%d) are not compatible.", nrow(x), length(y)))
  }

  # Create a list of indices of individuals for which some trait (or covariate) data is missing.
  miss = which(is.na(y))
  if ((hsreg.verbosity > 1) && (length(miss) > 0))
    cat("There are", length(miss), "subject(s) with missing trait data.\n")

  if (adjusted) 
  {
    miss.x = which(apply(is.na(x), 1, any))
    if ((hsreg.verbosity > 1) && (length(miss.x) > 0))
      cat("There are", length(miss.x), "subject(s) with missing covariate data.\n")
   
    miss = unique(c(miss, miss.x))
  }

  # Remove observations with missing data in trait or covariates.
  if (length(miss))
  {
    if (hsreg.verbosity > 0)
      cat("There are", length(miss), "subjects with missing trait or covariate data.\n")

    y = as.numeric(y[-miss])
    if (adjusted) 
      x = x[-miss, , drop = FALSE]
  }

  if (hsreg.verbosity > 0)
    cat(length(y), "subjects will be included in the analysis.\n")

  # Return the tidied-up y,x, geno, and removal list.
  return(list("y"=y, "x"=x, "missing"=miss))
}

##------------------------------------------------------------------------------------------------------
## Prepare the genotype data (by removing entries associated with missing data in y or x), and by
## converting it to a matrix of counts of "A alleles".
##------------------------------------------------------------------------------------------------------
hsreg.prepare.geno.data <- function(geno, geno.format, miss=NULL)
{
  # Remove any genotype data for which we have missing trait or covariate information.
  if (geno.format == "tped")
  {
    # TPED format: two columns per individual, with allele pairs.
    if (length(miss) > 0)
    {
      tmp = (miss-1)*2
      range.rm = c(rbind(tmp+1, tmp+2))
      geno = geno[,-range.rm]
    }

    # Get the allele counts.
    A = geno[,1]
    AA = (geno==A)
    ncols = ncol(AA)
    odds  = seq(1,ncols,2)
    evens = seq(2,ncols,2)
    AA.counts = AA[,odds] + AA[,evens]

    geno = AA.counts

  } else {

    # Impute format: 3 columns per individual, with posterior probabilities
    # for AA,Aa,aa genotypes.
    if ((ncol(geno) %% 3) != 0)
      stop("IMPUTE format genotype data should have 3 columns per individual.")

    if (length(miss) > 0)
    { 
      tmp = (miss-1)*3
      range.rm = c(rbind(tmp+1, tmp+2, tmp+3))
      geno = geno[,-range.rm]
    }

    # Convert the genotype data to counts of A alleles.
    n.cols = ncol(geno)
    AA.col = geno[,seq(1, n.cols, 3)]
    Aa.col = geno[,seq(2, n.cols, 3)]
    geno = 2*AA.col + Aa.col
  }

  # The algorithm actually uses 1 - geno.
  geno = 1 - geno
  
  return (geno)
}

##------------------------------------------------------------------------------------------------------
## hsreg.pmat.fun
## This is a helper function used to pre-calculate the P matrix and some associated values. These
## values all depend only on the trait and covariate data and therefore can be calculated once
## for a per-window analysis and then re-used for each window. This speeds up the p-value 
## calculations. 
##------------------------------------------------------------------------------------------------------
## y             = Vector of trait values, for a binary trait this must be 1 = case, 0 = control, no 
##                 missing values allowed.
## x             = Matrix of non-genetic covariates, EXCLUDING the intercept term.
##                 This can be NA if there is no covariate information, but if present, no missing 
##                 values are allowed.
## trait.type    = Character string defining the type of trait, with values of "gaussian" for a continuous 
##                 trait or "binomial" for binary trait.
##------------------------------------------------------------------------------------------------------
hsreg.pmat.fun <- function(y, x, trait.type)
{
  # Check the trait type specified.
  trait.int = hsreg.check.trait(trait.type)

  # Check for covariates, and whether the number of covariates matches the number of individuals.
  adjusted = !all(is.na(x))
  if (adjusted)
  {
    x <- as.matrix(x)
    if (nrow(x) != length(y)) 
      stop(gettextf("Dimensions of trait data (%d) and covariate data (%d) are not compatible.", nrow(x), length(y)))
  }

  n.subj <- length(y)

  if (adjusted)
  {
    ## Adjusted trait.
    P.mat.class = "A"    

    reg.out  = glm(y ~ x, family = trait.type)
    Ey = reg.out$fitted.values
    Vy = switch(trait.int, sum(reg.out$resid^2)/reg.out$df.resid, Ey*(1-Ey))
    newx = as.matrix( cbind(rep.int(1, n.subj), x) )

    if (trait.int == 1) # Adjusted Gaussian
    { 
      w = 1/Vy
      P.mat = diag(w, n.subj) - (w * newx) %*% tcrossprod( ginv(crossprod(newx)), newx)
    }
    else                # Adjusted Binomial
    {
      # Vy is a vector in this case.
      w = Vy
      wnewx = w * newx # Same as diag(w) %*% newx

      # P.mat = WW - WW %*% newx.adj %*% ginv(t(newx.adj) %*% WW %*% newx.adj ) %*% t(newx.adj) %*% WW
      # WW = diag(Vy)
      # crossprod(A,B) = t(A) %*% B
      # tcrossprod(A,B) = A %*% t(B)
      # So we have:
      # P.mat = W - W %*% X %*% ginv(t(W %*% X) %*% X) %*% t(W %*% X)
      # W = t(W) (since W is diagonal)
      # P.mat = W - W %*% X %*% ginv(t(X) %*% W %*% X) %*% t(X) %*% W

      P.mat = diag(w) - wnewx %*% tcrossprod( ginv( crossprod(wnewx, newx) ), wnewx )
    }
  }
  else
  {
    ## Unadjusted trait.
    P.mat.class = "U"    

    Ey = mean(y)
    Vy = var(y)

    if (trait.int == 1) # Unadjusted Gaussian
      w = 1/Vy
    else                # Unadjusted Binomial
      w = Vy

    if ((hsreg.ps.calc == 2) || (hsreg.ps.calc == 0)) # 0 = calibration time.
      P.mat = diag(w, n.subj) - w/n.subj
    else
      P.mat = NULL
  }

  newy = y - Ey
  if (trait.int == 1) newy = newy * 1/Vy

  return (list("pmat"=P.mat, "newy"=newy, "vy"=Vy, "pmat.class"=P.mat.class, "w"=w))
}

ps.calc.1 = function(geno, w)
{
  if (hsreg.timings) start = proc.time()
  nsubj = ncol(geno)
  nloci = nrow(geno)
  SSS = ((2/nloci) * crossprod(geno - 1)) + 2
  PS = (w * SSS)
  wcolsums = (w / nsubj) * colSums(SSS)
  for (i in 1:nsubj)
    PS[i,] = PS[i,] - wcolsums
  remove(SSS)
  remove(wcolsums)
  if (hsreg.timings) cat("PS (simple) time = ", (proc.time() - start)[3], "\n")
  return(PS)
}

ps.calc.2 = function(geno, p)
{
  # SSS = ((2/n.loci) * (t(geno) %*% (geno)) + 2
  # PS = P.mat %*% SSS

  if (hsreg.timings) start = proc.time()
  PS = (2/nrow(geno)) * tcrossprod(p, geno) %*% geno
  if (hsreg.timings) cat("PS (by geno) time =", (proc.time() - start)[3], "\n")

  # The second term is just the 2 * (row sums of P), but these are all zero.

  return(PS)
}

##------------------------------------------------------------------------------------------------------
## hsreg.pval.fun
##------------------------------------------------------------------------------------------------------
## y             = Vector of trait values, for a binary trait this must be 1 = case, 0 = control.
##                 No missing values allowed.
## geno          = The genotype data as prepared by hsreg.prepare.geno.
## x             = Matrix of non-genetic covariates, EXCLUDING the intercept term.
##                 This can be skipped if there is no covariate information. 
##                 No missing values allowed.
## trait.type    = Character string defining the type of trait, with values of "gaussian" for a continuous 
##                 trait or "binomial" for binary trait.
## approx.method = A number specifying the type of method used to approximate p-value, with values of
##                 2 = Gamma approximation (i.e., the 2-moment approximation), and
##                 3 = the 3-moment approximation
##                 Generally, the 2-moment approximation is considerably faster than the 3-moment.
## pmat.data       The P matrix and associated values as created by hsreg.prepare.trait.data
##------------------------------------------------------------------------------------------------------
hsreg.pval.fun <- function(y, geno, x=NA, trait.type="binomial", approx.method=2, pmat.data=NULL, in.tvc=NULL, in.ps=NULL)
{
  # Check the approximation method specified
  if ((approx.method != 2) && (approx.method != 3))
    stop(gettextf("Approx.method must be 2 (= 2-moment, gamma) or 3 (= 3-moment), but it is %s", approx.method))

  # Check the trait type specified.
  trait.int = hsreg.check.trait(trait.type)

  # If we have not been given the P matrix - calculate it.
  if (length(pmat.data) == 0)
    pmat.data = hsreg.pmat.fun(y,x,trait.type)
  P.mat = pmat.data$pmat
  P.mat.class = pmat.data$pmat.class
  newy = pmat.data$newy
  Vy = pmat.data$vy

  if (length(in.tvc) == 0) # We must calculate PS
  { 
    # Calculate VC
    # SSS = (2/n.loci) * (t(geno) %*% (geno)) + 2
     # Tvc = 1/2 * newy %*% SSS %*% newy
    if (hsreg.timings) start = proc.time()
    n.loci = nrow(geno)
    # snewy = sum(newy)
    tvc.tmp = tcrossprod(newy, geno)
#    Tvc = 1/n.loci * tcrossprod(newy, geno) %*% (geno %*% newy) +  2 * snewy * snewy  # CSCS newy = y - Ey, so is this always zero?
    Tvc = 1/n.loci * tcrossprod(tvc.tmp, tvc.tmp) # +  2 * snewy * snewy  # CSCS newy = y - Ey, so is this always zero?
    remove(tvc.tmp)
    if (hsreg.timings) cat("Tvc time =", (proc.time() - start)[3], "\n")

    if ((hsreg.ps.calc == 1) && (P.mat.class == "U"))
      PS = ps.calc.1(geno, pmat.data$w)
    else
      PS = ps.calc.2(geno, P.mat)
  }
  else
  {
    Tvc = in.tvc
    PS = in.ps
  }

  # Calculate p-value
  if (approx.method == 3)
  {
    ## 3-moment approximation
    if (length(in.tvc) == 0) # We were given the completed PS, so do nothing now.
    {
      w = pmat.data$w
      if (length(w) > 1)
      {
        r = diag(1/sqrt(w))
        PS = 1/2 * r %*% (PS %*% P.mat) %*% r
      }
      else
      {
        # w is just a scalar
        PS = 1/2 * 1/w * PS %*% P.mat
      }
    }
    # What we are calling PS here is actually PSP and is therefore symmetric.
    eg = eigen(PS, symmetric=T, only.values=T);
    evalue = eg$values; 
    ev = evalue[round(evalue,10)>0]
    c1 = sum(evalue);      
    c2 = sum(evalue^2);      
    c3 = sum(evalue^3);      
    dof = hprime = c2^3/c3^2
    pval = 1-pchisq( (Tvc-c1)*sqrt(hprime/c2)+hprime, dof)
  }
  else if (approx.method == 2) 
  {
    # Gamma (2-moment) approximation

    if (hsreg.timings) start = proc.time()
    trace.ps = sum(diag(PS))
    Evc = 1/2 * trace.ps
    #Itt = 1/2 * sum(diag(PS %*% PS)) 
    Itt = 1/2 * get.sum.diag.square(PS)

    # cat("Evc = ", Evc, " Itt = ", Itt, "\n");

    Vvc = Itt
    if (trait.int == 1) # Gaussian
    {
      # In the Gaussian case Vy is a single number.
      # Itp = 1/2 * sum(diag( Vy * PS ))/Vy^2
      # Ipp = 1/2 * sum(diag( P.mat * Vy ))/Vy^2
      Itp = 1/2 * trace.ps / Vy
      Ipp = 1/2 * sum(diag( P.mat )) / Vy
      Vvc = Itt - Itp^2/Ipp
    }

    if (hsreg.timings) cat("pval prep time =", (proc.time() - start)[3], "\n")

    # Get gamma.
    aa = Evc^2/Vvc
    bb = Evc/Vvc
    if (hsreg.timings) start = proc.time()
    pval = 1 - pgamma(Tvc, aa, bb)
    if (hsreg.timings) cat("pval time =", (proc.time() - start)[3], "\n")
  }

  result = list()
  result$tvc = Tvc
  result$pval = pval
  if (length(in.ps) == 0) result$ps = PS
  if (approx.method == 3) result$ev = ev ## JY 2/14/2011

  remove(PS)
  remove(P.mat)

  return (result)
}

## This function returns the sum of the diagonal elements in the square of a matrix.
## For a large matrix this is faster than calculating the entire square of the matrix
## and then summing the diagonals. 
get.sum.diag.square <- function(m)
{
  s = 0
  r = nrow(m)
  for (i in 1:r)
    s = s + drop(m[i,] %*% m[,i])
  return(s)
}

##------------------------------------------------------------------------------------------------------
## hsreg.vc.score.fun
##------------------------------------------------------------------------------------------------------
## y          = vector of trait values, for binary trait this must be 1 = case, 0 = control. 
##              No missing values allowed (see hsreg.prepare.trait.data).
## geno       = Genotype data as counts of A alleles for subject and locus.
## x          = Matrix of non-genetic covariates, EXCLUDING the intercept term.
##              This can be skipped if there is no covariate information. 
##              If present, no missing values allowed (see hsreg.prepare.trait.data).
## trait.type = Character string defining the type of trait, with values of "gaussian" for a continuous 
##              trait or "binomial" for binary trait.
## approx.method = A number specifying the type of method used to approximate p-value, with values of
##                 2 = Gamma approximation (i.e., the 2-moment approximation), and
##                 3 = the 3-moment approximation
##                 Generally, the 2-moment approximation is considerably faster than the 3-moment.
##------------------------------------------------------------------------------------------------------
hsreg.vc.score.fun <- function(y, geno, x=NA, trait.type="binomial", approx.method=2)
{
  geno = data.matrix(geno)

  # Remove any missing trait or covariate values, and the associated genotype data.
  prep.data = hsreg.prepare.trait.data(y,x,trait.type)
  y = prep.data$y
  x = prep.data$x
  miss = prep.data$missing

  # Prepare the genotype data as A allele counts.
  geno = hsreg.prepare.geno.data(geno, "impute", miss)

  # Calculate a p-value from the result.
  result = hsreg.pval.fun(y,geno,x,trait.type,approx.method)

  return(result)
}

hsreg.vc.score.fun.prep <- function(y, geno, x=NA, trait.type="binomial", approx.method=2)
{
  # Calculate a p-value from the result.
  result = hsreg.pval.fun(y,geno,x,trait.type,approx.method)

  return (result)
}

##------------------------------------------------------------------------------------------------------
## Read a config file (with lines of "name = value" format).
##------------------------------------------------------------------------------------------------------
## filename = name of config file to be read.
##------------------------------------------------------------------------------------------------------
hsreg.read.config.fun <- function(filename)
{
  lines = readLines(filename)

  config = list()
  config$use.covs = 1
  config$verbosity = 1
  config$save.psmat = 0
  config$out.dir = "."
  config$snp.start = 1
  config$geno.format = "impute"

  cat("\n")
  
  for (l in lines)
  {
    l = gsub("^ +", "", l)
    if (nchar(l) == 0) next;
    if (substr(l,1,1) == "#") next; # Comment

    cat("  ", l, "\n", sep="")
    
    splits = strsplit(l, "=")
    for (s in splits)
    {
      v = s[1] # First value is the name.
      lv = nchar(v)

      name = substring(l,1,lv)
      name = gsub("^ +", "", name)
      name = gsub(" +$", "", name)

      value = substring(l,lv+2)
      value = gsub("^ +", "", value)
      value = gsub(" +$", "", value)

      if (name == "win.size")
        config$win.size = as.numeric(value)
      else if (name == "win.start")
        config$win.start = as.numeric(value)
      else if (name == "snp.start")
        config$snp.start = as.numeric(value)
      else if (name == "snp.total")
        config$snp.total = as.numeric(value)
      else if (name == "geno.file")
        config$geno.file = value
      else if (name == "geno.format")
        config$geno.format = value
      else if (name == "geno.start.col")
        config$geno.start.col = as.numeric(value)
      else if (name == "trait.type")
        config$trait.type = value
      else if (name == "trait.file")
        config$trait.file = value
      else if (name == "use.covs")
        config$use.covs = as.numeric(value)
      else if (name == "approx.method")
        config$approx.method = value
      else if (name == "verbosity")
        config$verbosity = value
      else if (name == "save.psmat")
        config$save.psmat = as.numeric(value)
      else if (name == "save.smat") # Name was changed - this is for backwards compatability.
        config$save.psmat = as.numeric(value)
      else if (name == "print.timings")
        assign("hsreg.timings", as.numeric(value), .GlobalEnv)
      else if (name == "ps.calc")
        assign("hsreg.ps.calc", as.numeric(value), .GlobalEnv)
      else if (name == "run.id")
        config$run.id = value
      else if (name == "out.dir")
        config$out.dir = value
      else
        cat("Unexpected parameter name: ", name)
    }
  }

  return(config)
}

##------------------------------------------------------------------------------------------------------
## This function runs a windowed analysis.
##------------------------------------------------------------------------------------------------------
## filename = name of config file to be read.
## y        = trait vector
## x        = covariates
## miss     = rows in which some trait data is missing.
##------------------------------------------------------------------------------------------------------
hsreg.window.vc.score.fun <- function(config, y, x = NA)
{ 
  prep.data = hsreg.prepare.trait.data(y,x,config$trait.type)
  y = prep.data$y
  x = prep.data.x
  miss = prep.data$missing

  result = hsreg.vc.score.fun(config,y,x,miss)
}

hsreg.window.vc.score.fun.prep <- function(config, y, x = NA, miss = NULL)
{
  geno.file = config$geno.file
  geno.format = config$geno.format
  geno.start.col = config$geno.start.col 
  win.size = config$win.size
  snp.start = config$snp.start
  snp.total = config$snp.total
  trait.type = config$trait.type
  approx.method = config$approx.method
  out.dir = config$out.dir

  # Check the approximation method specified
  if ((approx.method != 2) && (approx.method != 3))
    stop(gettextf("Approx.method must be 2 (= 2-moment, gamma) or 3 (= 3-moment), but it is %s", approx.method))

  # Check the trait type specified.
  trait.int = hsreg.check.trait(trait.type)

  # Calculate number of windows.
  n.win = ceiling(snp.total / win.size)
  snp.max = snp.start + snp.total - 1

  results.file.name = paste("hsreg.window.results.txt", sep="")
  if (!is.null(config$run.id))
    results.file.name = paste("hsreg.window.results.", config$run.id, ".txt", sep="")
  results.file.name = paste(out.dir, "/", results.file.name, sep="")
  rfc = file(results.file.name, "w")
  cat(file=rfc, "SNPs Tvc Pval\n")
  
  # Open the genotypes file - we will read the required lines for each window as we go.
  gfc = file(geno.file, "r")

  # If we are not starting at the first SNP we must position the read pointer in the 
  # genotypes file to account for this.
  if (snp.start > 1)
  {
    skip = snp.start - 1
    chunk = 100
    while (skip > chunk) 
    {
      readLines(gfc, n=chunk)
      skip = skip - chunk
    }
    readLines(gfc, n=skip)
  }

  # Pre-calculate the P matrix, WW, newy and Vy
  pmat.data = hsreg.pmat.fun(y,x,trait.type)
  if (hsreg.verbosity > 9)
    cat("Prepared P matrix\n")

  # Now loop for each window in the genotypes file.
  start.time = proc.time()
  for (win in 1:n.win)
  {
    if (hsreg.verbosity > 0)
      cat("Window", win, "of", n.win, "\n")
    win.start = snp.start + (win.size * (win - 1))
    win.end = min(win.start + win.size - 1, snp.max)
    this.win.size = win.end - win.start + 1
    if (hsreg.verbosity > 1)
    {
      cat("- start SNP =", win.start, "\n")
      cat("- end SNP =", win.end, "\n")
      cat("- number of SNPs =", this.win.size, "\n")
    }

    # The genotype data is expected to be in IMPUTEv2 style. The first few columns
    # may contain non-genotype data. How many of these there are is indicated by 
    # geno.start.col.
    if (hsreg.timings) start = proc.time()
    geno = hsreg.read.geno(gfc, ignore=(geno.start.col-1), nsnps=this.win.size, geno.format=geno.format)
    if (hsreg.timings) cat("Read time =", (proc.time() - start)[3], "\n")
    
    # Format the genotype data as counts of A alleles.
    if (hsreg.timings) start = proc.time()
    geno = hsreg.prepare.geno.data(geno, config$geno.format, miss)
    if (hsreg.timings) cat("Geno prep time =", (proc.time() - start)[3], "\n")

    # Calculate the number of subjects (assumes IMPUTEv2 format)
    nsubj = ncol(geno)

    # Check the number of phenotypes specified.
    if (length(y) != nsubj) 
      stop(gettextf("Dimensions of trait data (%d) and genotype data (%d) are not compatible.", length(y), nsubj))

    # Check the end of the window.
    if (this.win.size != nrow(geno))
      stop(gettextf("Did not get the expected number of SNP rows from file: expected %d but got %d", this.win.size, nrow(geno)))

    # Calibration of PS calculation, for unadjusted case only.
    if ((win == 1) && (hsreg.ps.calc == 0) && (pmat.data$pmat.class == "U"))
    {
      start.time = proc.time()
      PS = ps.calc.1(geno, 1.0)
      ps.1.time = (proc.time() - start.time)[3]
      cat("PS 1 =", ps.1.time, "\n")

      start.time = proc.time()
      PS = ps.calc.2(geno, pmat.data$pmat)
      ps.2.time = (proc.time() - start.time)[3]
      cat("PS 2 =", ps.2.time, "\n")

      ps.calc = 1
      if (ps.2.time < ps.1.time) ps.calc = 2
      assign("hsreg.ps.calc", ps.calc, .GlobalEnv)
      cat("Will use PS calc", hsreg.ps.calc, "\n")
    }

    if (hsreg.timings) start = proc.time()
    result = hsreg.pval.fun(y=y, geno, x=x, trait.type=trait.type, approx.method=approx.method, pmat.data=pmat.data)
    if (hsreg.timings) cat("Pval time =", (proc.time() - start)[3], "\n")
    tvc = result$tvc
    pval = result$pval
    ps = result$ps
    cat(file=rfc, win, this.win.size, tvc, pval, "\n")

    if (hsreg.verbosity > 1)
    {
      cat("- Tvc =", tvc, "\n")
      cat("- p-value =", pval, "\n")
    }

    if (hsreg.timings) start = proc.time()
    ps = ps * this.win.size
    tvc = tvc * this.win.size

    if (win == 1)
    {
      cum.tvc = tvc
      cum.ps = ps
    }
    else
    {
      cum.tvc = cum.tvc + tvc
      cum.ps = cum.ps + ps
    }

    remove(ps) # Release memory
    if (hsreg.timings) cat("Cum time =", (proc.time() - start)[3], "\n")
  }

  close(rfc)
  close(gfc)

  end.time = proc.time()
  tm = end.time - start.time
  if (hsreg.verbosity > 0)
    cat("Average time per window =", round(tm[3] / n.win, 1), "secs\n")

  # Save the final PS (if requested). We save it before we divide by the total 
  # number of SNPs in this analysis making it easier to combine. 
  if (config$save.psmat)
  {
    ps.file.name = "hsreg.ps.txt"
    if (!is.null(config$run.id))
      ps.file.name = paste("hsreg.ps.",config$run.id,".txt",sep="")
    ps.file.name = paste(out.dir, "/", ps.file.name, sep="")
    
    # Write out the specifics of this analysis: SNP total, n.subj, window size
    cat(file=ps.file.name, snp.total, ncol(cum.ps), approx.method, cum.tvc, "\n")
    write.table(cum.ps, file=ps.file.name, quote=F, row.name=F, col.names=F, append=T)
  }

  # Correct for the total number of SNPs being used.
  cum.tvc = cum.tvc / snp.total
  cum.ps = cum.ps / snp.total

  # Get the overall p-value.
  result = hsreg.pval.fun(y=y, geno=NULL, x=x, trait.type=trait.type, approx.method=approx.method, in.tvc=cum.tvc, in.ps=cum.ps)

  # Save a file with the overall results.
  results.file.name = "hsreg.results.txt"
  if (!is.null(config$run.id))
    results.file.name = paste("hsreg.results.", config$run.id, ".txt", sep="")
  results.file.name = paste(out.dir, "/", results.file.name, sep="")
  final.result = c("SNPs"=snp.total, result$tvc, result$pval)
  write.table(file=results.file.name, t(final.result), quote=F)

  return(result)
}

##------------------------------------------------------------------------------------------------------
## This function reads genotype data from an impute format file.
##------------------------------------------------------------------------------------------------------
hsreg.read.geno = function(file, ignore=0, skip=0, nsnps=0, geno.format="dosage")
{
  if (is.character(file))
  {
    file = file(file, "rt")
    on.exit(close(file))
  }
  else
  {
    if (!inherits(file, "connection")) 
      stop("hsreg.read.geno: 'file' must be a character string or connection")
    if (!isOpen(file, "rt")) 
    {
      open(file, "rt")
      on.exit(close(file))
    }
  }

  # Read characters for tped data - else doubles.
  what = double(0)
  if (geno.format == "tped")
    what = character(0)

  hsreg.skip.geno(file, skip)

  s = 0
  while (s < nsnps)
  {
    # Discard initial columns as necessary.
    if (ignore > 1)
      scan(file, what=character(0), quiet=TRUE, nmax=ignore)

    # Read the rest of the row.
    row = scan(file, what=what, quiet=TRUE, nlines=1)

    # CSCS Specifying integer in the line above can save a small 
    # CSCS amount of time, but is only valid for straight 1/0 calls 
    # (not "dosage" data).
    # row = scan(file, what=integer(0), quiet=TRUE, nlines=1)

    if (s == 0)
      geno = matrix(row, 1, length(row))
    else
      geno = rbind(geno, row)
    s = s + 1
  }

  return(geno)
}

##----------------------------------------------------------------------------
## This function is called to skip to a start SNP in the genotype data file.
## The file parameter must be an open connection.
##----------------------------------------------------------------------------
hsreg.skip.geno = function(file, skip=0)
{
  if (!inherits(file, "connection")) 
    stop("hsreg.skip.geno: 'file' must be a connection")

  if (skip > 0)
  {
    chunk = 100
    while (skip > chunk) 
    {
      readLines(file, n=chunk)
      skip = skip - chunk
    }
    readLines(file, n=skip)
  }
}

##------------------------------------------------------------------------------------------------------
## This function runs an analysis based on the contents of a configuration file.
##------------------------------------------------------------------------------------------------------
## config.file = name of config file to be read, or a pre-populated config object.
## trait.data  = matrix of trait and covariate values (optional)
##
##               If trait.data is given it must be a matrix with the first column being identifiers, 
##               second column giving the trait values and subsequent columns (if any) giving covariate 
##               values. Missing values must be coded as NA.
##
##               If trait.data is not given, the config file must contain the name of a trait file 
##               suitable for reading with read.table. In the result the first column should be 
##               identifiers (not used), the second column must be the trait values, subsequent 
##               columns (if any) should be the covariate values. Missing values must be coded as 
##               NA. A binary trait should be coded as 0/1.               
##------------------------------------------------------------------------------------------------------
hsreg.analysis <- function(config, trait.data=NULL)
{
  hsreg.banner()

  if (is.null(config) || is.na(config))
    stop("you must specify a configuration file or object.")

  ptm = proc.time()

  # If the config specified is a character string, assume it is a file name.
  if (is.character(config))
  {
    cat("Reading config file:", config, "\n");
    if (!file.exists(config))
      stop(gettextf("config file specified does not exist: %s", config))
    config = hsreg.read.config.fun(config)
  }

  geno.file = config$geno.file
  geno.format = config$geno.format
  geno.start.col = config$geno.start.col 
  win.size = config$win.size
  snp.start = config$snp.start
  snp.total = config$snp.total
  trait.file = config$trait.file
  trait.type = config$trait.type
  approx.method = config$approx.method
  assign("hsreg.verbosity", config$verbosity, .GlobalEnv)
  
  windowed = !is.null(win.size)
  if (!windowed)
  {
    cat("\nAll SNPs analysis...\n")
  } else {
    cat("\nWindowed analysis...\n")
    cat("Window size =", win.size, "\n")
  }

  if (!is.null(config$run.id))
    cat("Run id =", config$run.id, "\n")

  cat("Calculation start SNP =", snp.start, "\n")
  cat("Calculation total SNPs =", snp.total, "\n")
  snp.max = snp.start + snp.total - 1
  cat("Calculation last SNP =", snp.max, "\n")
  cat("Genotypes file =", geno.file, "\n")
  cat("Geno start col =", geno.start.col, "\n")
  if (length(trait.file) > 0)
    cat("Traits file =", trait.file, "\n")
  cat("Trait type =", trait.type, "\n")
  cat("Approx method =", approx.method, "\n")

  # If trait data has not been given, read it from file. 
  if (length(trait.data) == 0)
    trait.data  = read.table(trait.file, header=T)

  # Extract the trait value column.
  # First column contains individual identifiers.
  y = trait.data[,2]

  # Extract the values of the explanatory variables (if any)
  x = NA
  num.cov = 0
  if (config$use.covs)
  {
    num.cov = ncol(trait.data) - 2
    if (num.cov > 0) 
      x  = trait.data[,c(3:(num.cov+2))]
  }

  prep.data = hsreg.prepare.trait.data(y,x,trait.type)
  y = prep.data$y
  x = prep.data$x
  miss = prep.data$miss

  if (!windowed)
  {
    # Load the genotype data.
    geno = hsreg.read.geno(geno.file, ignore=(geno.start.col-1), skip=snp.start-1, nsnps=snp.total, geno.format=geno.format)

    if (length(geno) == 0)
      stop("No genotype data read from file.")

    # Remove any genotype data for which we have missing trait or covariate information.
    geno = hsreg.prepare.geno.data(geno, config$geno.format, miss)

    # At this point the genotype data has been onverted to counts of A alleles per individual.

    # Run the all-in-one analysis.
    result = hsreg.vc.score.fun.prep(y=y, geno=geno, x=x, trait.type, approx.method=approx.method)

  } else {

    # Run the windowed analysis.
    result = hsreg.window.vc.score.fun.prep(config, y=y, x=x, miss=miss)
  }

  if (hsreg.verbosity > 0)
  {
    tm = proc.time() - ptm
    cat("Total time =", round(tm[3]), " secs\n")
  }

  return(result)
}


##------------------------------------------------------------------------------------------------------
## This function combines S matrices previously written out by a windowed analysis. Each ps file
## contains a header line giving the total number of SNPs in the run and the number of subjects used
## in the analysis (i.e. not including individuals with missing trait or covariate values).
##------------------------------------------------------------------------------------------------------
## ps.file.names = vector of ps files to be combined.
##------------------------------------------------------------------------------------------------------
hsreg.combine.ps.files = function(y, ps.file.names, x, trait.type="binomial", approx.method=2)
{
  if (length(ps.file.names) == 0)
    stop("no PS matrix file names specified")

  # Remove any missing trait or covariate values, and the associated genotype data.
  prep.data = hsreg.prepare.trait.data(y,x,trait.type)
  y = prep.data$y
  x = prep.data$x

  loci.total = 0
  for (ps.file.name in ps.file.names)
  {
    if (!file.exists(ps.file.name))
      stop(gettextf("file does not exist: %s", ps.file.name))

    # Read the analysis details
    ps.con = file(ps.file.name, "rt")
    details = scan(ps.con, n=4, quiet=T)
    n.loci = details[1]
    n.subj = details[2]
    this.approx.method = details[3]
    tvc = details[4]
    
    if (approx.method == 0)
      approx.method = this.approx.method
    else if (this.approx.method != approx.method)
      stop(gettextf("The requested approximation method, %s, and the method used to calculate the PS file, %s, do not match", approx.method, this.approx.method))

    # Read the PS matrix data.
    this.ps = matrix(scan(ps.con, n=n.subj*n.subj,quiet=T), nrow=n.subj, ncol=n.subj, byrow=T)
    close(ps.con)

    cat("Read ", ps.file.name, "\n")
    cat("n.loci =", n.loci, "n.subj =", n.subj, "Tvc =", tvc, "\n");

    if (loci.total == 0)
    {
      cum.tvc = tvc
      cum.ps = this.ps
    }
    else
    {
      cum.ps = cum.ps + this.ps
      cum.tvc = cum.tvc + tvc
    }

    # Keep track of the overall total number of SNPs.
    loci.total = loci.total + n.loci
  }

  # Correct the accummulated PS matrix for the total number of loci used.
  cum.tvc = cum.tvc / loci.total
  cum.ps = cum.ps / loci.total

  # Calculate the score statistic and p-value.
  result = hsreg.pval.fun(y,cum.ps,x,trait.type,approx.method, in.tvc=cum.tvc, in.ps=cum.ps)

  return (list("tvc"=result$tvc, "pval"=result$pval))
}
