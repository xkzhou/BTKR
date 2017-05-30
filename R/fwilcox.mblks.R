##----------------------------------------------------------------------------
## Wilcoxon rank-sum test for multiple experiments or blocks 2014/03/10 version
##----------------------------------------------------------------------------
fwilcox.mblks <- function(y, grp, blk){
  ## This function inplement the Wilcoxon rank-sum test for
  ## multiple experiments or blocks
  ## [Ref:Lehmann EL(1998), Nonparametrics: Statistical Methods Based on Ranks]
  ##browser()
  if(length(unique(grp))!=2)
    return("This function only compares difference in two groups")
  blk.uni <- unique(blk)
  Wi <- EWi <- VWi <- Ni <- rep(NA, length(blk.uni))
  for(i in 1:length(blk.uni)){
    id <- which(blk==blk.uni[i])
    r.sum <- tapply(rank(y[id], na.last="keep"),
                    grp[id], function(x) sum(x, na.rm=T))
    k <- which(r.sum==min(r.sum))[1]
    Wi[i] <- r.sum[k]
    n <- tapply(rank(y[id], na.last="keep"),
                grp[id], function(x) sum(!is.na(x)))
    Ni[i] <- sum(n)
    EWi[i] <- n[k]*(Ni[i]+1)/2
    VWi[i] <- prod(n)*(Ni[i]+1)/12
  }
  W <- sum(Wi/(Ni+1))
  EW <- sum(EWi/(Ni+1))
  VW <- sum(VWi/(Ni+1)^2)
  W.star <- (W-EW)/sqrt(VW)
  p.val <- 2*(1-pnorm(abs(W.star)))
  return(list(W=W.star, p.value=p.val))
}
