##-------------------------------------------------------------------
fpval.txt <- function(pval){
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
