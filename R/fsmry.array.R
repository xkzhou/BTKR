##-------------------------------------------------------------------------
fsmry.array <- function(dat.expr, dat.pheno, x.vars, base, dat.info){
  ## this function only summarize the high-throughput data by the dichotomized variables
  ## Note: the fold change below assumes that the metabolites levels
  ##       were log transformed
  out <- cbind(dat.info, dat.expr)

  ##browser()
  for (i in 1: length(x.vars)){
    library(qvalue)
    cat(paste(i, ".", x.vars[i], "\n", sep=""))
    x <- factor(dat.pheno[,x.vars[i]])
    ## y <-dat.expr[1,]
    out0 <- t(apply(dat.expr,1, function(y){
      tmp <- fsmry.by.grp(y, x, log.tr=F)
      tmp.avg <- tapply(y,x, mean)
      tmp.sd <- tapply(y,x, sd)
      tmp.n <- table(x)
      tmp.sigma <- sqrt(sum((tmp.n-1)*tmp.sd^2)/sum(tmp.n-1))
      tmp.efsize <- (tmp.avg[2]-tmp.avg[1])/tmp.sigma
      tmp.median <- tapply(y,x, median)
      tmp.fcavg <- base^(tmp.avg[2]-tmp.avg[1])
      tmp.fcmed <- base^(tmp.median[2]-tmp.median[1])
      tmp.pt <- as.numeric(attr(tmp, "detailed")[2,4])
      tmp.pw <- as.numeric(attr(tmp, "detailed")[2,6])
      tmp.pt.sig <- ifelse(tmp.pt<0.05,1,0)
      tmp.pt.sig[tmp.efsize<0&tmp.pt<0.05] <- -1
      tmp.pw.sig <- ifelse(tmp.pw<0.05,1,0)
      tmp.pw.sig[tmp.fcmed<1&tmp.pw<0.05] <- -1

      return(c(tmp.avg[1], tmp.sd[1], tmp.avg[2],tmp.sd[2],
               tmp.median, tmp.fcavg,tmp.efsize,
               tmp.pt,tmp.fcmed,tmp.pw, tmp.pt.sig, tmp.pw.sig))}))
    tmp.qt<-qvalue(as.numeric(out0[,9]))$qvalues
    tmp.qw<-qvalue(as.numeric(out0[,11]))$qvalues
    out0<-cbind(out0,tmp.qt,tmp.qw)
    dimnames(out0)[[2]] <-
      paste(c(rep(paste(gsub("\\.", "",x.vars[i]),levels(x),sep=""), each=2),
              paste(gsub("\\.", "",x.vars[i]),levels(x),sep=""),
              rep(x.vars[i],9)),
            c(rep(c(".avg", ".sd"),2), rep(".med",2),
              ".fcavg", ".efsize", ".pt", ".fcmed", ".pw",
              ".pt.sig", ".pw.sig",".qt",".qw"),
            sep="")

    out <- cbind(out, out0)
  }
  return(out)
}
