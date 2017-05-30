##-----------------------------------------------------------------------
##------------------------------------------------------------------------
## modified venn diaggram plotting function based on the function in ABarray
my.drawVennDiagram <- function (object, names, mar = rep(0.5, 4), cex = 1, ...)
{
  ## colors <- c("blue", "cyan", "orange")
  colors <- rep("black",3)
  if (!is(object, "VennCounts"))
    object <- vennCounts(object)
  nsets <- ncol(object) - 1
  if (nsets > 3)
    stop("Can't plot Venn diagram for more than 3 sets")
  if (missing(names))
    names <- colnames(object)[1:nsets]
  counts <- object[, "Counts"]
  totalCounts = sum(counts)
  theta <- 2 * pi * (1:360)/360
  xcentres <- list(0, c(-1, 1), c(-1, 1, 0))[[nsets]]
  ycentres <- list(0, c(0, 0), c(1/sqrt(5), 1/sqrt(5), -2/sqrt(5)))[[nsets]]
  r <- c(1.6, 1.6, 1.6)[nsets]
  xtext <- list(-1.2, c(-1.2, 1.2), c(-1.2, 1.2, 0))[[nsets]]
  ytext <- list(1.8, c(1.8, 1.8), c(2.4, 2.4, -3))[[nsets]]
  old.par <- par(mar = mar)
  on.exit(par(old.par))
  plot(x = 0, y = 0, type = "n", xlim = c(-4, 4), ylim = c(-4,
                                                           4), xlab = "", ylab = "", axes = FALSE, ...)
  for (circle in 1:nsets) {
    lines(xcentres[circle] + r * cos(theta), ycentres[circle] +
            r * sin(theta), col = colors[circle], lwd = 3)
    text(xtext[circle], ytext[circle], names[circle], cex = cex)
  }
  switch(nsets, {
    rect(-3, -2.5, 3, 2.5)
    text(2.3, -2.1, totalCounts, cex = cex)
    text(0, 0, counts[2], cex = cex)
  }, {
    rect(-3, -2.5, 3, 2.5)
    text(2.3, -2.1, totalCounts, cex = cex)
    text(1.5, 0.1, counts[2], cex = cex)
    text(-1.5, 0.1, counts[3], cex = cex)
    text(0, 0.1, counts[4], cex = cex)
  }, {
    rect(-3, -3.5, 3, 3.3)
    text(2.5, -3, totalCounts, cex = cex)
    text(0, -1.7, counts[2], cex = cex)
    text(1.5, 1, counts[3], cex = cex)
    text(0.75, -0.35, counts[4], cex = cex)
    text(-1.5, 1, counts[5], cex = cex)
    text(-0.75, -0.35, counts[6], cex = cex)
    text(0, 0.9, counts[7], cex = cex)
    text(0, 0, counts[8], cex = cex)
  })
  invisible()
}
my.doVennDiagram <- function (a, b, c = NULL, names, ...)
{
  require(limma)
  if (is.null(c)) {
    list.all <- union(a, b)
    list.mat <- matrix(0, nrow = length(list.all), ncol = 2)
    colnames(list.mat) <- c("list1", "list2")
    for (i in 1:length(list.all)) {
      list.mat[i, 1] <- list.all[i] %in% a
      list.mat[i, 2] <- list.all[i] %in% b
    }
    list.venn <- vennCounts(list.mat)
    my.drawVennDiagram(list.venn, names = names, cex = 1.5, ...)
    ab <- intersect(which(list.mat[, 1] == 1),
                    which(list.mat[, 2] == 1))
    list.ab <- vector(length = length(list.all))
    list.ab[ab] <- TRUE
    fileName <- "Venn_list_1+2.csv"
    if (!missing(names)) {
      fileName <- paste("Venn", names[1], names[2], ".csv",
                        sep = "_")
    }
    write.table(list.all[list.ab], file = fileName, sep = ",")
    print(paste("List information is written to file", fileName))
    invisible(list.all[list.ab])
  }
  else {
    list.all <- union(a, union(b, c))
    list.mat <- matrix(0, nrow = length(list.all), ncol = 3)
    colnames(list.mat) <- c("list1", "list2", "list3")
    for (i in 1:length(list.all)) {
      list.mat[i, 1] <- list.all[i] %in% a
      list.mat[i, 2] <- list.all[i] %in% b
      list.mat[i, 3] <- list.all[i] %in% c
    }
    list.venn <- vennCounts(list.mat)
    my.drawVennDiagram(list.venn, names = names, cex = 1.5, ...)
    ab <- intersect(which(list.mat[, 1] == 1), which(list.mat[,
                                                              2] == 1))
    ac <- intersect(which(list.mat[, 1] == 1), which(list.mat[,
                                                              3] == 1))
    bc <- intersect(which(list.mat[, 2] == 1), which(list.mat[,
                                                              3] == 1))
    abc <- intersect(which(list.mat[, 1] == 1), bc)
    vd <- matrix(0, nrow = length(list.all), ncol = 4)
    rownames(vd) <- list.all
    vd[, 1][ab] <- 1
    vd[, 2][ac] <- 1
    vd[, 3][bc] <- 1
    vd[, 4][abc] <- 1
    colnames(vd) <- c("list 1_2", "list 1_3", "list 2_3",
                      "list 1_2_3")
    if (!missing(names)) {
      colnames(vd) <- c(paste(names[1], names[2], sep = "+"),
                        paste(names[1], names[3], sep = "+"), paste(names[2],
                                                                    names[3], sep = "+"), paste(names[1], names[2],
                                                                                                names[3], sep = "+"))
    }
    fileName <- "Venn_list_1+2+3.csv"
    if (!missing(names)) {
      fileName <- paste("Venn", names[1], names[2], names[3],
                        ".csv", sep = "_")
    }
    write.table(vd, file = fileName, sep = ",", col.names = NA)
    print(paste("List information is written to file", fileName))
    invisible(vd)
  }
}
##----------------------------------------------------------------------
## a wrapper function of obtaining p-values for fixed effects estimate
## from lmer based on likelihood ratio test. This function is written by
## Christopher Moore
## http://blog.lib.umn.edu/moor0554/canoemoore/2010/09/lmer_p-values_lrt.html
## Example:
## library(lme4)
## library(SASmixed)
## lmer.out <- lmer(strength ~ Program * Time + (Time|Subj), data=Weights)
## p.values.lmer(lmer.out)

p.values.lmer <- function(x) {
  summary.model <- summary(x)
  data.lmer <- data.frame(model.matrix(x))
  names(data.lmer) <- names(fixef(x))
  names(data.lmer) <- gsub(pattern=":", x=names(data.lmer),
                           replacement=".", fixed=T)
  names(data.lmer) <- ifelse(names(data.lmer)=="(Intercept)",
                             "Intercept", names(data.lmer))
  string.call <- strsplit(x=as.character(x@call), split=" + (", fixed=T)
  var.dep <- unlist(strsplit(x=unlist(string.call)[2], " ~ ", fixed=T))[1]
  vars.fixef <- names(data.lmer)
  formula.ranef <- paste("+ (", string.call[[2]][-1], sep="")
  formula.ranef <- paste(formula.ranef, collapse=" ")
  formula.full <- as.formula(paste(var.dep, "~ -1 +",
                                   paste(vars.fixef, collapse=" + "),
                                   formula.ranef))
  data.ranef <- data.frame(x@frame[,
                                   which(names(x@frame) %in% names(ranef(x)))])
  names(data.ranef) <- names(ranef(x))
  data.lmer <- data.frame(x@frame[, 1], data.lmer, data.ranef)
  names(data.lmer)[1] <- var.dep
  out.full <- lmer(formula.full, data=data.lmer, REML=F)
  p.value.LRT <- vector(length=length(vars.fixef))
  for(i in 1:length(vars.fixef)) {
    formula.reduced <- as.formula(paste(var.dep, "~ -1 +",
                                        paste(vars.fixef[-i],
                                              collapse=" + "),
                                        formula.ranef))
    out.reduced <- lmer(formula.reduced, data=data.lmer, REML=F)
    print(paste("Reduced by:", vars.fixef[i]))
    print(out.LRT <- data.frame(anova(out.full, out.reduced)))
    p.value.LRT[i] <- round(out.LRT[2, 7], 3)
  }
  summary.model@coefs <- cbind(summary.model@coefs, p.value.LRT)
  summary.model@methTitle <- c("\n", summary.model@methTitle,
                               "\n(p-values from comparing nested models fit by maximum likelihood)")
  print(summary.model)
}
##---------------------------------------------------------------
## The following code is udsed for obtaing mcmcpvalues for mcmc samples from lmer object. It was provided by Dr. Douglas Bates:
##http://rwiki.sciviews.org/doku.php?id=guides:lmer-tests

mcmcpvalue <- function(samp)
{
  ## elementary version that creates an empirical p-value for the
  ## hypothesis that the columns of samp have mean zero versus a
  ## general multivariate distribution with elliptical contours.

  ## differences from the mean standardized by the observed
  ## variance-covariance factor
  std <- backsolve(chol(var(samp)),
                   cbind(0, t(samp)) - colMeans(samp),
                   transpose = TRUE)
  sqdist <- colSums(std * std)
  sum(sqdist[-1] > sqdist[1])/nrow(samp)
}
##------------------------------------------------------------------------
##The following function is used to merge multiple metabolon metabolomic data
fmet.merge <- function(datnames=dat.names, cols.ann, by.col="COMP_ID",
                       names.ids,
                       names.mets,
                       mets.start){
  ##browser()
  comps <- nc.mets <- NULL
  nc.ann <- length(cols.ann)
  n.dat <- length(datnames)

  for (i in 1:n.dat){
    x <- get(datnames[i])
    nc.mets <- c(nc.mets, dim(x)[2]-mets.start[i]+1)
    cols.na <- cols.ann[is.na(match(cols.ann, names(x)))]
    if(length(cols.na)>0){
      x.add <- data.frame(matrix(NA,nrow=dim(x)[1], nc=length(cols.na)))
      dimnames(x.add)[[2]] <- cols.na
      x <- data.frame(x.add,x)
    }
    comps <- rbind(comps,x[,cols.ann])
  }
  comp.ann <- apply(comps, 2, function(x){
    tapply(x, as.vector(comps[, by.col]), function(y)
      ifelse(any(!is.na(y)),y[!is.na(y)][1], NA))})

  dat.mets <- data.frame(matrix(NA,nrow=dim(comp.ann)[1], ncol=sum(nc.mets)))
  dimnames(dat.mets)[[2]] <- names.mets
  nc.end <- cumsum(nc.mets)
  nc.start <- c(1, nc.end[-length(nc.end)]+1)
  dat.id <- data.frame(matrix(0,nrow=dim(comp.ann)[1], ncol=n.dat))
  dimnames(dat.id)[[2]] <- names.ids
  for( i in 1:n.dat){
    x <- get(datnames[i])
    id <- match(x[,by.col], as.vector(comp.ann[,by.col]))
    dat.id[id,i] <- 1
    dat.mets[id, nc.start[i]:nc.end[i]] <- x[,mets.start[i]:dim(x)[2]]
  }
  dat.all <- data.frame(dat.id, comp.ann, dat.mets)
  return(dat.all)
}
##-------------------------------------------------------------
##Pwer and sample size calculation based on log-rank test
##

fpower.logrank <- function(n,p, hr, sig.level, power){
  ## this function calculate the power and sample size based on log-rank test for survival data
  ## n=number of events
  ## p=proportion of subjects in the experimental group
  ## hr=hazard ratio

  ## Power for a log-rank test is the same as power for a Cox model.
  ## For a 2-arm study
  ##          d = (qnorm(.975) + qnorm(.85))^2 / (.2 *.8 * coef^2)

  ##          d = # deaths required

  ##           .975 = two sided alpha=.05
  ##           .85   = power of .85
  ##            .2, .8 = proportion in each group

  ## coef = Cox model coef you want power against.
  ##   So for a 50% difference in hazard rates (or median survival times)
  ##   coef = log(1.5).
  ##  For a 50% change I get d=341. Now for the hard part of sample size
  ##  in a survival study: how many people to you need to enroll and how
  ##  long will you need to follow them, to observe 341 total deaths? This
  ##  second step is usually a mix of prior knowlege, enrollment expectations,
  ##  and wild speculation. In the words of a prior director of research at
  ##  U of Rochester:"At the commencement of a study the incidence of the
  ##  disease in question will drop by one half, and will not return to its
  ##  former levels until the study ends or the principle investigator retires,
  ##  whichever comes first." L.L.
  ##     -- Terry Therneau

  flr <- function(n,p, hr, sig.level,power){
    n-(qnorm(sig.level/2)+qnorm(power))^2/(p*(1-p)*log(hr)*log(hr))
  }


}


##  Do power calculation of Freedman for Cox PH as set forth on page 733
##  of 5th edition of Rosner, "Fundamentals of Biostatistics"
##
##            k = ratio of (# subjects in exposed group) to
##                         (# subjects in control)
##            t = max time of follow-up
##           RR = (hazard rate for exposed group) /
##                (hazard rate for control group)
##    lambda[j] = Rosner's lambda_(j-1)
##     delta[j] = Rosner's delta_(j-1)
##

survivalUtil.power <- function(n.1, n.2, t, RR, lambda, delta, alpha=0.05) {
  pf <- survivalUtil.probFail(t, RR, lambda, delta)
  k <- n.1/n.2
  m <- n.1*pf$p.E + n.2*pf$p.C
  pnorm(sqrt(k*m)*abs(RR - 1)/(k*RR + 1) - qnorm(1 - alpha/2))
}


##
##  Do sample-size calculation of Freedman for Cox PH as set forth
##  on page 735
##  of 5th edition of Rosner, "Fundamentals of Biostatistics"
##
##            k = ratio of (# subjects in exposed group) to
##                         (# subjects in control)
##            t = max time of follow-up
##           RR = (hazard rate for exposed group) /
##                (hazard rate for control group)
##    lambda[j] = Rosner's lambda_(j-1)
##     delta[j] = Rosner's delta_(j-1)
##
##  Value: n1 - number of people needed in exposed group
##         n2 - number of people needed in control group
##
survivalUtil.sampSize <- function(k, t, RR, lambda, delta,
                                  alpha=0.05, beta=.2) {
  m <- ((1/k)*(((k*RR + 1)/(RR - 1))^2)*
          ((qnorm(1 - (alpha/2)) + qnorm(1 - beta))^2))
  pf <- survivalUtil.probFail(t, RR, lambda, delta)
  n.1 <- (m*k)/((k*pf$p.E) + pf$p.C)
  n.2 <- m/((k*pf$p.E) + pf$p.C)
  list(n.1=n.1, n.2=n.2)
}


##
##  Compute probabilties of failure used by survivalUtil.power() and
##  survivalUtil.sampSize()
##
survivalUtil.probFail <- function(t, RR, lambda, delta) {
  A <- cumprod(1 - lambda)[1:t]
  B <- cumprod(1 - (RR*lambda))[1:t]
  C <- cumprod(1 - delta)[1:t]
  D <- lambda[2:(t + 1)]*A*C
  E <- RR*lambda[2:(t + 1)]*B*C
  list(p.C=sum(D), p.E=sum(E))
}

##------------------------------------------------------------------
## Sample size calculation for diagnostic test based on Pepe MS (2003),
## The StattisticalEvaluation of Medical Tests for Classification and Prediction,
## Oxford University Press, New York.
fn.dtest <- function(d0.t, d1.t, d0.f, d1.f, TPFb, FPFb, alpha, beta, paired){
  ## This function calculate the sample size needed for comparing two diagnostic
  ## tests with gold standards available
  if(paired){
    TPPF <- (1+d1.t)*TPFb-1
    FPPF <- (1+d1.f)*FPFb-1
    if(TPPF<0) TPPF <- 0
    if(FPPF<0) FPPF <- 0
    ##number of cases
    n.d <- ((qnorm(sqrt(1-alpha))+qnorm(sqrt(1-beta)))/log(d1.t/d0.t))^2*
      (((d1.t+1)*TPFb-2*TPPF)/(d1.t*TPFb^2))
    ## number of controls needed
    n.nd <- ((qnorm(sqrt(1-alpha))+qnorm(sqrt(1-beta)))/log(d1.f/d0.f))^2*
      (((d1.f+1)*FPFb-2*FPPF)/(d1.f*FPFb^2))
  }
  else{
    n.d <- ##number of cases
      n.d <- ((qnorm(sqrt(1-alpha))+qnorm(sqrt(1-beta)))/log(d1.t/d0.t))^2*
      ((d1.t+1-2*d1.t*TPFb)/(d1.t*TPFb))
    ## number of controls needed
    n.nd <- ((qnorm(sqrt(1-alpha))+qnorm(sqrt(1-beta)))/log(d1.f/d0.f))^2*
      ((d1.f+1-2*d1.f*FPFb)/(d1.f*FPFb))
  }
  return(c(n.d=n.d, n.nd=n.nd))
}

##--------------------------------------------------------------
## Bayesian stopping rules for phase II trials with binary outcomes
##---------------------------------------------------------------
fp.tox <- function(r0=0.1, a0, b0, n.max=40){
  ## this function calculate the posterior probability that the toxicity
  ## rate exceeds the prespecified level r0 given various observations for
  ## sample size up to n.max, prior belief of the success/failure
  ## rate is specified using a beta distribution beta(a0, b0)
  smry.prior <- list(mean=a0/(a0+b0),
                     W.90=c("5%"=qbeta(0.05,a0,b0),
                            "95%"=qbeta(0.95,a0,b0)))
  smry.post <- vector(length=n.max, mode="list")
  for(i in 1:n.max){
    x <- 0:i
    smry.post[[i]] <- 1-pbeta(r0, a0+x, b0+i-x)
    names(smry.post[[i]]) <- as.character(0:i)
  }
  return(list(smry.prior=smry.prior, smry.post=smry.post))
}
