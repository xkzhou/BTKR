##-------------------------------------------------------------------
fpval.txt <- function(x, d=3, p=TRUE){
  cp <- 10^(-1*d)
  id <- which(x<cp)
  out<- format(round(x,d), digits=d, nsmall=d)
  if(p) out <- paste0("P=", out)
  if(length(id)>0){
    out[id] <- paste0("<", cp)
    if(p) out[id] <- paste0("P<", cp)
  }
  return(out)
}

fpval2.txt <- function(pval, digits=3){
  ## this function summarize the p-value
  pval <- as.numeric(pval)
  if(pval<0.001) pval.txt <- "P<0.001"
  else if(pval<0.01) pval.txt <- paste("P=", round(pval,3),sep="")
  else if(pval>0.045 &pval<0.05)pval.txt <- "P<0.05"
  else pval.txt <- paste("P=", round(pval,2),sep="")

  return(pval.txt)
}

fpval.sbl <- function(pval){
  ## this function summarize the p-value using symbols
  pval <- as.numeric(pval)
  if(pval<0.001) pval.txt <- "***"
  else if(pval<0.01) pval.txt <- "**"
  else if(pval<0.05) pval.txt <- "*"
  else pval.txt <- ""

  return(pval.txt)
}
fauc.txt <- function(ci.auc, digits=2){
  ## ci.auc is obtained using ci(auc(outcome, pred)) using library(pROC)
  auc <- as.numeric(ci.auc)
  out <-paste("AUC=", round(auc[2], digits=digits),
              ", 95% CI: ", round(auc[1], digits=digits), "-",
              round(auc[3], digits=digits), sep="")
  return(out)
}

##-------------------------------------------
fcphuni.stat <- function(y=out.coxph){
  ## extract important outcomes
  n <- y$n
  nevent <-y$nevent

  y.smry <- summary(y)
  yname <- all.vars(y$formula)[2]
  xname <- all.vars(y$formula)[-(1:2)]
  xlevs <- y$xlevels
  y.stat <- cbind(y.smry$conf.int[,-2, drop=F], y.smry$coef[, 5, drop=F])
  p.lr <- y.smry$sctest["pvalue"]
  if(is.null(xlevs)){
    y.stat <- data.frame(levels="unit", y.stat,row.names = NULL)
    xref <- NULL
  } else if(y$contrasts[[1]]=="contr.poly"){
    y.stat <- data.frame(levels=row.names(y.stat), y.stat,row.names = NULL)
    xref <- NULL
  } else {
    y.stat <- data.frame(Variable=xlevs[[1]][-1], y.stat,row.names = NULL)
    xref <- xlevs[[1]][1]
  }
  colnames(y.stat) <- c("Level", "HR", "Lower95", "Upper95", "P.Wald")
  out <-list(stats=y.stat, p.lr= p.lr,
             yname=yname, xname=xname, xref=xref)
  return(out)
}

fcphuni.tbl <- function(x=out.fcphuni.stat){
  xlen <-length(x)
  out <- list()
  for (i in 1:xlen){
    xstat <- x[[i]]$stats
    rx.tbl <- cbind("X(levs)"=as.vector(xstat[,1]),
                    data.frame(lapply(xstat[,2:4], function(x)
                      format(round(x,2), digits=3, nsmall=2))),
                    P.Wald= fpval.txt(xstat[,5]),
                    P.LR=rep("", dim(xstat)[1]))
    if(i==1)
      r1.tbl <- c(paste0("**", x[[i]]$xname, "**"), rep("", 4),
                  as.vector(fpval.txt(x[[i]]$p.lr)))
    else if(i==2)
      r1.tbl <- c(paste0("*", gsub("x\\.", "",x[[i]]$xname), "*"),
                  rep("", 4), as.vector(fpval.txt(x[[i]]$p.lr)))
    r1.tbl <- data.frame(t(r1.tbl))
    names(r1.tbl) <- names(rx.tbl)
    if(is.null(x[[i]]$xref))
      out[[i]] <- rbind(r1.tbl, rx.tbl)
    else{
      r2.tbl <- c(x[[i]]$xref, "reference", rep("", 4))
      r2.tbl <- data.frame(t(r2.tbl))
      names(r2.tbl) <- names(rx.tbl)
      out[[i]] <- rbind(r1.tbl, r2.tbl, rx.tbl)
    }
  }
  out <-do.call(rbind, out)
  return(out)
}
