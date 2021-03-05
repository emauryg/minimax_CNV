##------------------------------------------------------------------------------
## ped.to.impute
##------------------------------------------------------------------------------
## Convert PLINK .ped style genotypes to IMPUTEv2 style.
## We assume that PLINK style data has 2 columns per marker, one row 
## per individual. Each pair of entries per marker per individual indicates 
## the alleles present in that individual.
## We generate IMPUTEv2 style data with each row representing one marker and
## each set of three consecutive columns representing the posterior probability
## of the genotype (AA,Aa,aa). Since we started with specific genotype calls in
## the PLINK format data we fill the impute matrix in with 1s and 0s only. 
## This function works on an in-memory matrix and will be useful for small
## to medium sized data sets. 
##------------------------------------------------------------------------------
snp.to.impute<-function(SNP)
{
  nloci = ncol(SNP)
  nsubj = nrow(SNP)
  impute = matrix(0, ncol=(3 * nsubj), nrow=nloci)
  for(row in 1:nloci) 
  {
    for (i in 1:nsubj)
    {
      genotype = SNP[row,i]
      col = 3*(i-1) + 1
      impute[row,col+genotype] = 1
    }
  }
  return(impute)
}



ped.to.impute<-function(ped)
{
  ncols = ncol(ped)
  if (ncols %% 2) 
    stop("Ped genotype data must have an even number of columns.") 

  ## Get one allele at each locus.
  A = vector(length = ncols)
  loci = seq(1, ncols, 2)
  for(j in loci) 
  {
    alleles = unique(c(ped[,j], ped[,j+1]))

    ## We can only cope with two alleles per marker
    if (length(alleles) > 2) 
      stop("More than two alleles for column ", j)

    ## We just need to choose A for this marker.
    A[j] = alleles[1]
  }

  nrows = nrow(ped)
  impute = matrix(0, ncol=(3 * nrows), nrow=ncols/2)
  for(j in loci) 
  {
    row = (j+1)/2
    for (i in 1:nrows)
    {
      genotype = 0 # 0=aa, 1=Aa, 2=AA
      if (ped[i,j] != ped[i,j+1]) 
        genotype = 1 # Must be Aa (or aA) 
      else if (ped[i,j] != A[j]) # A from above.
        genotype = 2

      col = 3*(i-1) + 1
      impute[row,col+genotype] = 1
    }
  }

  return(impute)
}

##------------------------------------------------------------------------------
## haplo.to.impute
##------------------------------------------------------------------------------
## Convert haplo.score style data to IMPUTEv2 style. We are expecting as input just the
## genotype data - no preceding columns with family id, individual id etc.
haplo.to.impute <- function(haplo)
{
  return (ped.to.impute(haplo))
}

