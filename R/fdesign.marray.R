##--------------------------------------------------------------------------
##Microarray study design
##--------------------------------------------------------------------------
fdesign.marray <- function(n=20, a1=0.5, a2=0.5, m=1000, m1=50, fdr=0.05,
                           d=1.0, r1){
  ## This function calculates the number of truly differentially expressed
  ## genes to be observed (r1) at the 5% fdr cutoff, given the effect size is
  ## delta, assuming we have n arrays, na1 in group 1 nad na2 in group2,
  ## and the array consists of a total of m genes and m1 are true
  ## differentially expressed with effect size of delta when a two sided test
  ## is used, and r1 is the number of true rejections. This function is based on formula 10 in (Sin-Ho Jung,
  ## Bioinformatics, Vol. 21 (no. 14), 2005)


  alpha <- (r1*fdr)/((m-m1)*(1-fdr))
  beta <- 1-r1/m1
  return(n-1-(qnorm(1-alpha/2)+qnorm(1-beta))^2/(a1*a2*d^2))

  ## may use the following to find the answer for one of the variables,
  ## for example r1:
  ## uniroot(fdesign.marray,interval=c(0.1,50), n=20, a1=0.5, a2=0.5, m=1000,
  ##  m1=50, fdr=0.05, d=2.0)$root
}
