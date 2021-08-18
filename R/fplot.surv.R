#' Plot Kaplan-Meier Curves with number at risk and number censored at bottom
#'
#' @param smdl.yt survival time
#' @param smdl.ye event
#' @param smdl.x a covariate
#' @param xnames  group names of the covariate x
#' @param xlab label of x-axis
#' @param ylab label of y-axis
#' @param dat data
#' @param t.at location of the x-axis tick marks
#' @param mar margin parameters
#' @param mgp mgp parameters
#' @param nr.line lines from the x-axis for the number at risk texts
#' @param nr.xadj adjustment along x-axis for the NR(NC) notation as a fraction of the length of the second tick mark from 0
#' @param nr.yadj adjustment along y-axis for the NR(NC) notation
#' @param loc.legend location of the legend
#' @param loc.p location of the p-value
#'
#' @return survfit output, median survival and coxph output if there is a covariate x
#'
#' @examples
#' set.seed(16)
#' dat.work <- data.frame(StudyID=1:40,
#'                        Death = sample(c(0,1), 40, replace=T, p=c(0.2, 0.8)),
#'                        Mths.OS = 10*rexp(40, rate=0.5),
#'                        Group = factor(rep(c("grp1", "grp2"), c(15,25))))
#' out <- fplot.surv(smdl.yt = "Mths.OS",
#'                   smdl.ye = "Death",
#'                   xlab = "Months",
#'                   ylab = "Overall Survival (%)",
#'                   dat = dat.work,
#'                   mar=c(4.0,3.5,0.5,0.7),
#'                   mgp = c(2.0, 0.8, 0),
#'                   nr.line = 3,
#'                   nr.xadj = 0.4,
#'                   nr.yadj = 0.8)
#'
#' out <- fplot.surv(smdl.yt = "Mths.OS",
#'                   smdl.ye = "Death",
#'                   smdl.x = "Group",
#'                   xnames = c("Group1", "Group2"),
#'                   xlab = "Months",
#'                   ylab = "Overall Survival (%)",
#'                   dat=dat.work,
#'                   t.at = NULL,
#'                   lty=NULL,
#'                   mar = c(4.0,3.5,0.5,0.7),
#'                   mgp = c(2.0, 0.8, 0),
#'                   nr.line = 3,
#'                   nr.xadj = 0.4,
#'                   nr.yadj = 0.8,
#'                   loc.legend = "bottomleft",
#'                   loc.p = "topright")
#'
fplot.surv <- function(smdl.yt = "Mths.OS",
                       smdl.ye = "OS",
                       smdl.x = "",
                       xnames="",
                       xlab = "Months",
                       ylab = "Overall Survival (%)",
                       dat=dat.work,
                       t.at = NULL,
                       lty=1,
                       mar=c(4.0,3.5,0.5,0.5),
                       mgp=c(2.0, 0.8, 0),
                       nr.line=3,
                       nr.xadj=0.75,
                       nr.yadj=0.8,
                       loc.legend="bottomleft",
                       loc.p="topright"){
  ##require(survival)

  smdl <- paste0("Surv(", smdl.yt, ", ", smdl.ye, ") ~ ", ifelse(smdl.x=="",1,smdl.x))
  if(is.null(t.at)) t.at <- pretty(dat[, smdl.yt])
  yt.max <- max(dat[, smdl.yt])
  t.sfit <- t.at
  tat.max <-max(t.at)
  if(tat.max > yt.max) t.sfit[length(t.sfit)] <- yt.max
  if(is.null(lty)) lty <- 1:length(xnames)

  out.sfit <- survfit(formula(smdl), data=dat)
  smry.sfit <- summary(out.sfit, time=t.sfit)

  mar[1] <- mar[1]+(0.2+nr.yadj)*(length(xnames)-1)
  par(mar=mar, mgp=mgp)
  plot(out.sfit, mark.time = T, lty=lty, xlim=range(t.at),
       lwd=2,las=1,font=2,
       xlab=xlab, ylab=ylab, axes=F,font.lab=2, cex.lab=0.8)
  axis(1, at=t.at, cex.axis=0.7, font.axis=2)
  axis(2, las=1, cex.axis=0.7, font.axis=2)
  box(lwd=2)

  med.surv <-NULL

  if(smdl.x==""){
    n.risk <- smry.sfit$n.risk
    n.event <- smry.sfit$n.event
    n.censor <- smry.sfit$n.censor
    if(tat.max > yt.max) {
      n <-length(t.at)
      n.risk[n] <- n.risk[n-1]-n.event[n]-n.censor[n]
    }
    n.censor <- cumsum(n.censor)
    mtext("NR (NC)", side=1, line=nr.line, at=-nr.xadj*t.at[3], outer=F, font=2, cex=0.7, adj=0)
    mtext(paste0(n.risk, "(", n.censor, ")"),
          side=1, line=nr.line, at=t.at, font=2, cex=0.7)
    med.surv <- rbind(med.surv, unlist(quantile(out.sfit, p=0.5)))
    colnames(med.surv) <- c("Med.Surv", "Lower.95CI", "Upper.95CI")
  } else {
    x.levels <- levels(smry.sfit$strata)
    out.coxph <- coxph(formula(smdl), data=dat)
    pval.txt <- BTKR::fpval.txt(summary(out.coxph)$sctest["pvalue"])

    legend(loc.legend, legend = xnames, lty=lty, bty="n", cex=0.8)
    legend(loc.p, legend=pval.txt, cex=0.8, bty="n")

    mtext("NR (NC)", side=1, line=nr.line-nr.yadj, at=-nr.xadj*t.at[3], outer=F, font=2, cex=0.7, adj=0)
    for(i in 1:length(x.levels)){
      id <- which(smry.sfit$strata==x.levels[i])
      n.risk <- smry.sfit$n.risk[id]
      n.event <- smry.sfit$n.event[id]
      n.censor <- smry.sfit$n.censor[id]
      tat.id <- t.at[1:length(id)]
      yt.max <- max(smry.sfit$time[id])
      if(max(tat.id) > yt.max) {
        n <-length(t.at)
        n.risk[n] <- n.risk[n-1]-n.event[n]-n.censor[n]
      }
      n.censor <- cumsum(n.censor)
      mtext(xnames[i], side=1, line=nr.line+nr.yadj*(i-1), at=-nr.xadj*t.at[3], outer=F, font=2, cex=0.7, adj=0)
      mtext(paste0(n.risk, "(", n.censor, ")"),
            side=1, line=nr.line+nr.yadj*(i-1), at=tat.id, font=2, cex=0.7)
      med.surv <- rbind(med.surv, unlist(quantile(out.sfit[i], p=0.5)))
    }
    colnames(med.surv) <- c("Med.Surv", "Lower.95CI", "Upper.95CI")
    row.names(med.surv) <- x.levels
  }
  if(length(xnames)>1) out <-list(sfit=out.sfit, medSurv=med.surv, coxph=out.coxph)
  else out <- list(sfit=out.sfit, medSurv=med.surv)
  return(out)
}
