#' Summarize a continuous variable by a categorical variable
#'
#' @param y A continuous variable.
#' @param grp A categorical group variable.
#' @param log.tr An indicator of whether to log transform y
#' @param method Text string for multiple comparisons method if there are more than 2 categories, such as "Tukey", "Dunnet", or specifying the contrast
#' @param method.padj Text string for specifying method for multiple adjustment of p-values, e.g., "bonferroni", "fdr", etc.
#' @param paired An indicator of whether paired test should be used
#' @param IQR An indicator of whether to summarize y using IQR (default) or range.
#' @return Summary of difference in y by grp using parametric and non-parametric methods.
#' @examples
#' fsmry.by.grp(y=c(rnorm(10), rnorm(10, mean=0.5)),
#'              grp=factor(rep(c("grp1", "grp2"), each=10)))
fsmry.by.grp <- function(y, grp,
                         log.tr=F,
                         method="Tukey",
                         method.padj=NULL,
                         paired=F,
                         IQR=T){
  ## note that measurements in different groups are assumed independent here
  ## random effects models may need to be used
  ## require library(multcomp)
  ## source("~/packages/my_R_code.R")
  library(multcomp)
  if (log.tr) y <- log(y)
  if(!is.factor(grp)) grp <- factor(grp)
  n <- table(grp)
  ##browser()
  n.cmplt <- tapply(y, grp, function(x) sum(!is.na(x)))
  mean.sd <- tapply(y, grp, function(x)
    paste(round(mean(x, na.rm=T),2), " +/- ",
          round(sd(x, na.rm=T),2),sep=""))
  if(IQR)
    median.IQR <- tapply(y, grp, function(x)
      paste(round(median(x, na.rm=T),2), " (",
            round(quantile(x, p=0.25, na.rm=T),2),",",
            round(quantile(x, p=0.75, na.rm=T),2),")",sep=""))
  else
    median.range <- tapply(y, grp, function(x)
      paste(round(median(x, na.rm=T),2), " (",
            round(min(x, na.rm=T),2),",",
            round(max(x, na.rm=T),2),")",sep=""))


  if(length(levels(grp))>2){
    pval.aov <- anova(lm(y~grp))[1,"Pr(>F)"]
    ##browser()
    ##summary(glht(lm(y~grp), linfct=mcp(grp=method)))
    multcmp.cname <- NULL
    if(is.null(method.padj)) method.padj <- "holm"
    if(length(method)>1){
      if(method.padj!="none"){
        multcmp.out <-
          data.frame(summary(glht(lm(y~grp),
                                  linfct=mcp(grp=method)),
                             test=adjusted(method.padj))$test[3:6])
        multcmp.cname <-
          c(multcmp.cname,
            paste("pm.p.", method.padj, sep=""))

        multcmp.out <-
          cbind(multcmp.out,
                data.frame(summary(glht(lm(y~grp),
                                        linfct=mcp(grp=method)),
                                   test=adjusted("none"))$test[3:6]))
      }
      else
        multcmp.out <- data.frame(summary(glht(lm(y~grp),
                                               linfct=mcp(grp=method)),
                                          test=adjusted("none"))$test[3:6])
      multcmp.cname <- c(multcmp.cname,"pm.p")
    }
    else{
      multcmp.out <-
        data.frame(summary(glht(lm(y~grp),
                                linfct=mcp(grp=method)))$test[3:6])
      multcmp.cname <- c(multcmp.cname, paste("pm.p.",method, sep=""))
    }
    names.comp <- rownames(multcmp.out)

    pval.KW <- kruskal.test(y~grp)$p.value
    p.wilcox <- rep(NA, length(names.comp))
    if(length(method)==1) p.t <- rep(NA, length(names.comp))
    ## browser()
    for(i in 1:length(names.comp)){
      grp.names <- strsplit(names.comp[i], split=" - ")[[1]]
      y1 <- y[grp==grp.names[1]]
      y2 <- y[grp==grp.names[2]]
      p.wilcox[i] <- wilcox.test(y1,y2, exact=F)$p.value
      if(length(method)==1){
        if(length(y1)==1&length(y2)==1)
          p.t[i] <- 1
        else if(length(y1)==1|length(y2)==1)
          p.t[i] <- t.test(y1,y2, var.equal=T)$p.value
        else
          p.t[i] <- t.test(y1,y2, var.equal=F)$p.value
      }
    }
    p.wilcox.adj <- p.adjust(p.wilcox,method=method.padj)
    if(length(method)==1){
      if(method.padj=="none"){
        multcmp.out <- data.frame(multcmp.out,p.t, p.wilcox)[,4:6]
        multcmp.cname <- c(multcmp.cname, "pm.p", "wilcoxon.p")
      }
      else{
        multcmp.out <- data.frame(multcmp.out,p.t,p.wilcox.adj, p.wilcox)[,4:7]
        multcmp.cname <-
          c(multcmp.cname, "pm.p",paste("wilcoxon.p.",method.padj, sep=""),
            "wilcoxon.p")
      }
    }
    else{
      if(method.padj=="none"){
        multcmp.out <- data.frame(multcmp.out, p.wilcox)[,4:5]
        multcmp.cname <- c(multcmp.cname, "wilcoxon.p")
      }
      else{
        multcmp.out <-
          data.frame(multcmp.out, p.wilcox.adj, p.wilcox)[,c(4,8:10)]
        multcmp.cname <-
          c(multcmp.cname, paste("wilcoxon.p.",method.padj,sep=""),
            "wilcoxon.p")
      }
    }
    dimnames(multcmp.out)[[2]] <- multcmp.cname
    pvalue.pairwise.comp <- data.frame(apply(multcmp.out,2, function(x)
      ifelse(x<=0.001,"<0.001",round(x,3))))
    ##    dimnames(pvalue.pairwise.comp)[[1]] <- names.comp
    ##    dimnames(pvalue.pairwise.comp)[[2]] <- c("pm.padj", "Wilcoxon.padj")
    if(IQR){
      out <- data.frame(n=as.numeric(n),
                        n.complete=n.cmplt,
                        mean.sd=mean.sd,
                        pval.anova=c(rep("",length(n)-1),
                                     ifelse(pval.aov<0.001, "<0.001", round(pval.aov,3))),
                        median.IQR=median.IQR,
                        pval.KW=c(rep("",length(n)-1),
                                  ifelse(pval.KW<0.001, "<0.001", round(pval.KW,3))))
      out2 <- data.frame(n=as.numeric(n),
                         n.complete=n.cmplt,
                         mean.sd=mean.sd,
                         pval.anova=c(rep("",length(n)-1),pval.aov),
                         median.IQR=median.IQR,
                         pval.KW=c(rep("",length(n)-1),pval.KW),
                         stringsAsFactors=F)
    }
    else{
      out <- data.frame(n=as.numeric(n),
                        n.complete=n.cmplt,
                        mean.sd=mean.sd,
                        pval.anova=c(rep("",length(n)-1),
                                     ifelse(pval.aov<0.001, "<0.001", round(pval.aov,3))),
                        median.range=median.range,
                        pval.KW=c(rep("",length(n)-1),
                                  ifelse(pval.KW<0.001, "<0.001", round(pval.KW,3))))
      out2 <- data.frame(n=as.numeric(n), n.complete=n.cmplt,
                         mean.sd=mean.sd,
                         pval.anova=c(rep("",length(n)-1),pval.aov),
                         median.range=median.range,
                         pval.KW=c(rep("",length(n)-1),pval.KW),
                         stringsAsFactors=F)
    }
    attr(out, "detailed") <- list(summary=out2, pval.pairwise.comp=multcmp.out)
    return(list(summary=out,pval.pairwise.comp=pvalue.pairwise.comp))
  }
  else{
    ##browser()
    grp.levels <- levels(grp)
    pval.t <- t.test(y[grp==grp.levels[1]], y[grp==grp.levels[2]],
                     paired=paired, na.action=na.omit)$p.val
    pval.w <- wilcox.test(y[grp==grp.levels[1]], y[grp==grp.levels[2]],
                          paired=paired, na.action=na.omit)$p.val
    if(IQR){
      out <- data.frame(n=as.numeric(n), n.complete=n.cmplt,
                        mean.sd=mean.sd,
                        pval.t=c(rep("",length(n)-1),
                                 ifelse(pval.t<0.001, "<0.001", round(pval.t,3))),
                        median.IQR=median.IQR,
                        pval.W=c(rep("",length(n)-1),
                                 ifelse(pval.w<0.001, "<0.001", round(pval.w,3))))
      out2 <- data.frame(n=as.numeric(n), n.complete=n.cmplt,
                         mean.sd=mean.sd,
                         pval.t=c(rep("",length(n)-1),pval.t),
                         median.IQR=median.IQR,
                         pval.W=c(rep("",length(n)-1),pval.w), stringsAsFactors=F)
    }
    else{
      out <- data.frame(n=as.numeric(n), n.complete=n.cmplt,
                        mean.sd=mean.sd,
                        pval.t=c(rep("",length(n)-1),
                                 ifelse(pval.t<0.001, "<0.001", round(pval.t,3))),
                        median.range=median.range,
                        pval.W=c(rep("",length(n)-1),
                                 ifelse(pval.w<0.001, "<0.001", round(pval.w,3))))
      out2 <- data.frame(n=as.numeric(n), n.complete=n.cmplt,
                         mean.sd=mean.sd,
                         pval.t=c(rep("",length(n)-1),pval.t),
                         median.range=median.range,
                         pval.W=c(rep("",length(n)-1),pval.w), stringsAsFactors=F)
    }
    attr(out,"detailed") <- out2
    return(out)
  }
}
