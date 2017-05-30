fsmry2.by.grp <- function(y=dat.shs[,"gender"],
                          grp=dat.shs[,"case_control"],
                          cmp.method=c("fisher","chisq")){
  ## This function summarize information for categorical variables
  ##browser()
  n.y <- table(y)
  n.missing <- sum(is.na(y))
  if (n.missing>0)
    n.mis.tab <- table(ifelse(is.na(y),"Yes", "No"))

  y.by.grp <- table(y, grp)

  if(cmp.method=="fisher"){
    test.out <- fisher.test(y.by.grp)
    if (n.missing>0){
      mis.by.grp <- table(ifelse(is.na(y),"Yes", "No"), grp)
      test2.out <- fisher.test(mis.by.grp)
    }
  }
  if(cmp.method=="chisq"){
    test.out <- chisq.test(y.by.grp)
    if(n.missing>0){
      mis.by.grp <- table(ifelse(is.na(y),"Yes", "No"), grp)
      test2.out <- chisq.test(mis.by.grp)
    }
  }

  out.y <- data.frame(n=as.numeric(n.y),
                      apply(y.by.grp,2, function(x)
                        paste(x," (",round(100*x/sum(x),2),"%)", sep="")),
                      p.value=c(rep("",length(n.y)-1),
                                ifelse(test.out$p.value<0.001,"<0.001",
                                       round(test.out$p.value,3))))
  row.names(out.y) <- row.names(y.by.grp)
  attr(out.y,"detailed") <- test.out
  if(n.missing>0){
    out.missing <- data.frame(n=as.numeric(n.mis.tab),
                              apply(mis.by.grp,2, function(x)
                                paste(x," (",round(100*x/sum(x),2),"%)",
                                      sep="")),
                              p.value=c(rep("",length(n.mis.tab)-1),
                                        ifelse(test2.out$p.value<0.001,"<0.001",
                                               round(test2.out$p.value,3))))
    row.names(out.missing) <- row.names(mis.by.grp)
    attr(out.missing,"detailed") <- test2.out
    return(list(summary=out.y, missingness=out.missing))
  }
  else
    return(out.y)
}
