#' Summarize a categorical variable by a categorical variable
#'
#' @param y A categorical variable.
#' @param grp A categorical group variable.
#' @param  An indicator of whether to log transform y
#' @param cmp.method Text string for comparisons method, "fisher" or, "chisq"
#' @param prop.by.row Logic value of whether to calculate the proportions by row, default is FALSE
#' @return Summary of difference in distribution of y by grp using Fisher's exact test or chi-squared test.
#' @examples
#' fsmry2.by.grp(y = factor(rep(c("Female", "Male", "Female", "Male"), c(2,8,6,4))),
#'               grp = factor(rep(c("grp1", "grp2"), each=10)),
#'               cmp.method = "fisher")
#' fsmry2.by.grp(y = factor(rep(c("Female", "Male", "Female", "Male"), c(2,8,6,4))),
#'               grp = factor(rep(c("grp1", "grp2"), each=10)),
#'               cmp.method = "fisher", prop.by.row=T)
fsmry2.by.grp <- function(y, grp, cmp.method=c("fisher","chisq"), prop.by.row=FALSE){
  ## This function summarize information for categorical variables
  ## browser()
  y <- y[!is.na(grp)]
  grp <- grp[!is.na(grp)]
  n.y <- table(y)
  n.missing <- sum(is.na(y))
  if (n.missing>0)
    n.mis.tab <- table(ifelse(is.na(y),"Yes", "No"))

  y.by.grp <- table(y, grp)

  if(cmp.method=="fisher"){
    test.out <- try(fisher.test(y.by.grp), TRUE)
    if(class(test.out)=="try-error")
      test.out <- fisher.test(y.by.grp, simulate.p.value = T)
    if (n.missing>0){
      mis.by.grp <- table(ifelse(is.na(y),"Yes", "No"), grp)
      test2.out <- try(fisher.test(mis.by.grp), TRUE)
      if(class(test2.out)=="try-error")
        test2.out <- fisher.test(mis.by.grp, simulate.p.value = T)
    }
  }
  if(cmp.method=="chisq"){
    test.out <- chisq.test(y.by.grp)
    if(n.missing>0){
      mis.by.grp <- table(ifelse(is.na(y),"Yes", "No"), grp)
      test2.out <- chisq.test(mis.by.grp)
    }
  }

  if(prop.by.row){
    y.by.grp.tbl <- t(apply(y.by.grp, 1, function(x)
      paste(x," (",round(100*x/sum(x),2),"%)", sep="")))
    } else {
    y.by.grp.tbl <- apply(y.by.grp, 2, function(x)
      paste(x," (",round(100*x/sum(x),2),"%)", sep=""))
  }
  dimnames(y.by.grp.tbl) <- dimnames(y.by.grp)

  out.y <- data.frame(n=as.integer(n.y),
                      y.by.grp.tbl,
                      p.value=c(rep("",length(n.y)-1),
                                ifelse(test.out$p.value<0.001,"<0.001",
                                       round(test.out$p.value,3))))
  attr(out.y,"detailed") <- test.out

  if(n.missing>0){
    if(prop.by.row){
      mis.by.grp.tbl <- t(apply(mis.by.grp, 1, function(x)
        paste(x," (",round(100*x/sum(x),2),"%)", sep="")))
    } else {
      mis.by.grp.tbl <- apply(mis.by.grp, 2, function(x)
        paste(x," (",round(100*x/sum(x),2),"%)", sep=""))
    }
    dimnames(mis.by.grp.tbl) <- dimnames(mis.by.grp)

    out.missing <- data.frame(
      n=as.integer(n.mis.tab),
      mis.by.grp.tbl,
      p.value=c(rep("",length(n.mis.tab)-1),
                ifelse(test2.out$p.value<0.001,"<0.001", round(test2.out$p.value,3))))
    attr(out.missing,"detailed") <- test2.out
    return(list(summary=out.y, missingness=out.missing))
  }
  else
    return(out.y)
}
