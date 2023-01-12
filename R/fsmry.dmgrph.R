#' Generating summary statistics of data in a data frame by a categorical variable
#'
#' @param dat A data frame.
#' @param vars A vector of variable names to be analyzed.
#' @param vars.cat A vector of indicators of whether the variable is categorical or not
#' @param vars.chisq A vector of indicators specifying for which variables the chi-square test should be used
#' @param by A categorial variable name the vars will be summarized by.
#' @param all An indicator of whether to provide the overall summary statistics of the data
#' @param markdown An indicator of whether to generate the summary table by bolding the variable names using markdown language.
#' @param IQR An indicator of whether to summarize the continuous using IQR (default) or range.
#' @param prop.by.row Logic value of whether to calculate the proportions by row, default is FALSE
#'
#' @return Summary statistics of each variable in a data frame and by grp using parametric and non-parametric methods.
#'
#' @examples
#' set.seed(16)
#' dat.work <- data.frame(ht = c(rnorm(10, mean=1.72, sd=0.1), rnorm(10, mean=1.65, sd=0.1)),
#'                        wt = c(rnorm(10, mean=70, sd=10), rnorm(10, mean=60, sd=10)),
#'                        sex = factor(rep(c("Female", "Male", "Female", "Male"), c(2,8,6,4))),
#'                        group = factor(rep(c("grp1", "grp2"), each=10)))
#'
#' fsmry.dmgrph(dat = dat.work,
#'              vars = c("ht", "wt", "sex"),
#'              vars.cat = c(0, 0, 1),
#'              by =  "group",
#'              prop.by.row=F)
#'
fsmry.dmgrph <- function(dat=dat.work,
                         vars=vars,
                         vars.cat=vars.cat,
                         vars.chisq=rep(0,length(vars)),
                         by="BMI.cat",
                         all=T,
                         markdown=T,
                         IQR=T,
                         prop.by.row=FALSE){
  ## This function produce the demographics data for all study subjects
  ## and/or by a categorical variable
  ## Example: out.byCohort <- fsmry.dmgrph(dat=dat.work,vars=vars,
  ##                                       vars.cat=vars.cat,
  ##                                       by="cohort", all=T)
  ## write.table(out.byCohort[[1]], "DemograpicTable__by_cohort_20121212.txt",
  ##             sep="\t",row.name=F, col.name=T, quote=T)

  ##browser()

  if (class(dat)[1] == "tbl_df")
    dat <- data.frame(dat)

  id.na <- which(is.na(dat[,by]))
  if(length(id.na)>0) dat <- dat[-id.na,]
  n <- dim(dat)[1]
  out <- NULL
  if(all){
    out.all <- vector(length=length(vars), mode="list")
    names(out.all) <- vars
  }
  if(!is.null(by)){
    out.by.grp <- vector(length=length(vars), mode="list")
    names(out.by.grp) <- paste(vars," ~ ", by ,sep="")
    dat[,by] <- factor(dat[,by])
    level.by <- levels(dat[,by])
    n.level <- length(level.by)
    n.by.cat <- table(dat[,by])
  }
  for (i in 1:length(vars)){
    cat(paste(i, ". ", vars[i],"\n", sep=""))
    tmp <- NULL
    ##browser()
    if(vars.cat[i]){ # for a categorical variable
      if(all){
        if(sum(is.na(dat[,vars[i]]))==0)
          out.all[[i]] <- cbind(n=table(dat[,vars[i]]),
                                prop=prop.table(table(dat[,vars[i]])))
        else{
          n.missing <- sum(is.na(dat[,vars[i]]))
          out.all[[i]] <-
            cbind(n=c(table(dat[,vars[i]]),
                      missing=n.missing),
                  prop=c(prop.table(table(dat[,vars[i]])),
                         missing=n.missing/n))
        }
        tmp <-
          rbind(c(ifelse(markdown,
                         paste("**", vars[i], ", n (%)**",sep=""),
                         paste(vars[i], ", n (%)",sep="")), ""),
                cbind(paste("   ", row.names(out.all[[i]])),
                      paste(out.all[[i]][,"n"], "(",
                            round(out.all[[i]][,"prop"]*100,2),"%)",
                            sep="")))
      }
      if(!is.null(by)){
        if(!vars.chisq[i])
          out.by.grp[[i]] <- fsmry2.by.grp(y=dat[,vars[i]],
                                           grp=dat[,by],
                                           cmp.method="fisher",
                                           prop.by.row=prop.by.row)
        else{
          ##browser()
          out.by.grp[[i]] <- fsmry2.by.grp(y=dat[,vars[i]],
                                           grp=dat[,by],
                                           cmp.method="chisq",
                                           prop.by.row=prop.by.row)
        }
        if(is.null(tmp)){
          if(sum(is.na(dat[,vars[i]]))==0)
            tmp <-
              rbind(c(ifelse(markdown,
                             paste("**",vars[i], ", n (%)**",sep=""),
                             paste(vars[i], ", n (%)",sep="")),
                      rep("", n.level+1)),
                    as.matrix(cbind(paste("   ",row.names(out.by.grp[[i]]), sep=""),
                                    as.matrix(out.by.grp[[i]][,-1]))))
          else
            tmp <-
              rbind(c(ifelse(markdown,
                             paste("**",vars[i], ", n (%)**",sep=""),
                             paste(vars[i], ", n (%)",sep="")),
                      rep("", n.level+1)),
                    as.matrix(cbind(c(row.names(out.by.grp[[i]]), "   missing"),
                                    rbind(as.matrix(out.by.grp[[i]][[1]][,-1]),
                                          as.vector(out.by.grp[[i]][[2]][2,-1])))))
        }
        else {
          if(sum(is.na(dat[,vars[i]]))==0)
            tmp <- rbind(c(tmp[1,],rep("", n.level+1)),
                         as.matrix(cbind(tmp[-1,],
                                         as.matrix(out.by.grp[[i]][,-1]))))
          else
            tmp <-
              rbind(c(tmp[1,], rep("", n.level+1)),
                    as.matrix(cbind(tmp[-1,],
                                    rbind(as.matrix(out.by.grp[[i]][[1]][,-1]),
                                          as.vector(out.by.grp[[i]][[2]][2,-1])))))
        }
      }
      out <- rbind(out,tmp)
    }
    else{ #for continuous variable
      if(all){
        if(IQR)
          out.all[[i]] <- c(mean=mean(dat[,vars[i]], na.rm=T),
                            sd=sd(dat[,vars[i]], na.rm=T),
                            median=median(dat[,vars[i]], na.rm=T),
                            iqr.L=quantile(dat[,vars[i]], p=0.25, na.rm=T),
                            iqr.U=quantile(dat[,vars[i]], p=0.75, na.rm=T))
        else
          out.all[[i]] <- c(mean=mean(dat[,vars[i]], na.rm=T),
                            sd=sd(dat[,vars[i]], na.rm=T),
                            median=median(dat[,vars[i]], na.rm=T),
                            min=min(dat[,vars[i]], na.rm=T),
                            max=max(dat[,vars[i]], na.rm=T))
        tmp <- rbind(c(ifelse(markdown,
                              paste("**",vars[i],"**",sep=""),
                              vars[i]), ""),
                     c("   Mean+/-sd",
                       paste(round(out.all[[i]][1],2),"+/-",
                             round(out.all[[i]][2],2),sep="")),
                     c(ifelse(IQR, "   Median (IQR)",
                              "   Median (range)"),
                       paste(round(out.all[[i]][3],2)," (",
                             round(out.all[[i]][4],2),", ",
                             round(out.all[[i]][5],2), ")",sep="")))
        if(sum(is.na(dat[,vars[i]]))>0){
          n.missing <- sum(is.na(dat[,vars[i]]))
          out.all[[i]] <- c(out.all[[i]],
                            n.missing=n.missing,
                            prop.missing=n.missing/n)
          tmp <- rbind(tmp,
                       c("   missing, n (%)",
                         paste(out.all[[i]]["n.missing"], " (",
                               round(100*out.all[[i]]["prop.missing"],2), "%)",
                               sep="")))
        }
      }
      if(!is.null(by)){
        out.by.grp[[i]] <- fsmry.by.grp(y=dat[,vars[i]],
                                        grp=dat[,by],
                                        log.tr=F, IQR=IQR)
        if(length(level.by)==2) obg <- out.by.grp[[i]]
        else obg <- out.by.grp[[i]][[1]]
        if(is.null(tmp)){
          tmp <- rbind(c(ifelse(markdown,
                                paste("**",vars[i],"**",sep=""),
                                vars[i]),
                         rep("",n.level+1)),
                       c("   Mean+/-sd", as.vector(obg[,3]),
                         as.vector(obg[n.level,4])),
                       c(ifelse(IQR, "   Median (IQR)",
                                "   Median (range)"),
                         as.vector(obg[,5]),
                         as.vector(obg[n.level,6])))
          if(sum(is.na(dat[,vars[i]]))>0){
            n.missing <- obg[,1]-obg[,2]
            prop.missing <- n.missing/obg[,1]
            tmp <- rbind(tmp,
                         c("   missing, n (%)",
                           paste(as.vector(n.missing), " (",
                                 as.vector(round(100*prop.missing,2)),
                                 "%)", sep="")))
          }
        }
        else{
          tmp1 <- rbind(c(tmp[1,],rep("",n.level+1)),
                        c(tmp[2,], as.vector(obg[,3]),
                          as.vector(obg[n.level,4])),
                        c(tmp[3,], as.vector(obg[,5]),
                          as.vector(obg[n.level,6])))
          if(dim(tmp)[1]>3) {
            if(!vars.chisq[i])
            #if(any(table(ifelse(is.na(dat[,vars[i]]),"Y","N"), dat[,by])<=5))
              mis.by.grp <- fsmry2.by.grp(y=ifelse(is.na(dat[,vars[i]]),"Y","N"),
                                          grp=dat[,by],
                                          cmp.method="fisher")
            else
              mis.by.grp <- fsmry2.by.grp(y=ifelse(is.na(dat[,vars[i]]),"Y","N"),
                                          grp=dat[,by],
                                          cmp.method="chisq")
            tmp1 <- rbind(tmp1,
                          c(tmp[4,],as.vector(unlist(mis.by.grp[2,-1]))))
          }
          tmp <- tmp1
        }
      }
      out <- rbind(out,tmp)
    }
  }
  ##browser()
  out <- data.frame(out)
  if(all){
    dimnames(out)[[2]][1:2] <- c("Variables",
                                 paste("All(n=",n,")", sep=""))
    if(!is.null(by))
      dimnames(out)[[2]][3:(2+n.level+1)] <-
        c(paste(level.by, rep("(n=",n.level),
                n.by.cat, rep(")", n.level),sep=""),
          "p.value")
  }
  else {
    dimnames(out)[[2]][1] <- "Variables"
    if(!is.null(by))
      dimnames(out)[[2]][2:(1+n.level+1)] <-
        c(paste(level.by, rep("(n=",n.level),
                n.by.cat, rep(")", n.level),sep=""),
          "p.value")
  }
  if(!is.null(by))
    return(list(demographic=out,summary.all=out.all,
                summary.by.grp=out.by.grp))
  else
    return(list(demographic=out,summary.all=out.all))
}
