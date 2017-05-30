##---------------------------------------------------------------
## A wrapper function for summarizing the univariate association
## between x and y and generating the corresponding graphs
##---------------------------------------------------------------
fsmry.graph <- function(y, x,type=c("scatter","bxp","errbar","bar","bar2"),
                        fname=NULL,subgrp="",stat.txt=NULL,loc.stat=NULL,
                        y.fnlab=NULL, y.plab,
                        x.fnlab=NULL, x.plab,xnames=NULL,
                        withline=F,geomean=F,
                        width=3, height=3,ylim=NULL,
                        mar=NULL, mgp=NULL,
                        cex.bar=NULL,
                        cex=1,
                        ...){
  ## This function is used to generate commonly used plots
  ##browser()
  if(!is.null(fname)){
    if(fname=="sys"){
      fdate <- gsub("-", "", Sys.Date())
      file.name <- paste("Fig/fig_",type, "_", gsub("\\.","_",y.fnlab),
                         "_",gsub("\\.","_",x.fnlab),"_",subgrp,"_",
                         fdate,".wmf", sep="")
    }
    else
      file.name <- fname
    win.metafile(file=file.name,width=width, height=height)
  }

  if(type=="scatter"){
    if(is.null(mar))mar <- c(3.5,3.5,0.5,0.5)
    if(is.null(mgp)) mgp <- c(2.3,0.8,0)
    par(mar=mar, mgp=mgp)
    ylim <- range(y, na.rm=T)
    xlim <- range(x, na.rm=T)
    plot(x,y, xlab=x.plab, ylab=y.plab,pch=20, cex=cex, las=1,...)
    if(withline)
      abline(coef(line(x,y)),...)
    if(!is.null(stat.txt)){
      if(is.null(loc.stat))
        text(xlim[1], ylim[2],label=stat.txt, cex=cex, pos=4)
      else
        text(loc.stat[1], loc.stat[2],label=stat.txt, cex=cex)
    }
  }
  if(type=="bxp"){
    if(is.null(mar))mar <- c(2.5,3.5,0.5,0.5)
    if(is.null(mgp)) mgp <- c(2.2,0.8,0)
    par(mar=mar, mgp=mgp)
    if(is.null(ylim))
      y.lim <- range(pretty(c(floor(min(y, na.rm=T)),
                              ceiling(max(y, na.rm=T)))))
    else y.lim <- ylim
    if(is.null(xnames)) xnames <- levels(x)
    boxplot(y~x, boxwex=0.5, xlab="", ylab=y.plab, names=xnames, ylim=y.lim,
            las=1, cex.axis=1,...)
    if(!is.null(stat.txt)){
      if(is.null(loc.stat))
        text(x=0.9, y=y.lim[2]-0.01*(y.lim[2]-y.lim[1]),label=stat.txt, cex=cex)
      else
        text(loc.stat[1], loc.stat[2],label=stat.txt, cex=cex)
    }
  }
  if(type=="errbar"){
    ##library(Hmisc)
    if(is.null(mar))mar <- c(2.5,3.5,0.5,0.5)
    if(is.null(mgp)) mgp <- c(2.2,0.8,0)
    par(mar=mar, mgp=mgp)
    if(geomean){
      y.mean <- tapply(y, x,function(z) exp(mean(log(z), na.rm=T)))
      y.sd <- tapply(y,x,function(z) exp(sd(log(z), na.rm=T)))
    }
    else{
      y.mean <- tapply(y, x,function(z) mean(z, na.rm=T))
      y.sd <- tapply(y,x,function(z) sd(z, na.rm=T))
    }

    y.lim <- range(0, ceiling(y.mean+1.5*y.sd))
    if(is.null(xnames)) xnames <- levels(x)
    mp <- barplot(y.mean,
                  space=c(0.2,rep(0.4, length(xnames)-1))/0.6,
                  width=0.6,
                  xlim=c(0,length(xnames)),col="black", names="",
                  ylim=y.lim, xlab="", ylab=y.plab, las=1)
    segments(x0=mp, x1=mp, y0=y.mean, y1=y.mean+y.sd,...)
    segments(x0=mp-0.025, x1=mp+0.025, y0=y.mean+y.sd, y1=y.mean+y.sd,...)
    ##errbar(x=mp, y=y.mean, yplus=y.mean+y.sd, yminus=y.mean-y.sd,
    ##       pch=".", add=T)
    axis(1, at=mp, labels=xnames, tick=0,...)
    axis(2, las=1,...)
    box(...)
    if(is.null(loc.stat))
      text(x=mp[1]-0.1, y=y.lim[2]-0.05*y.lim[2],label=stat.txt, cex=cex)
    else
      text(loc.stat[1], loc.stat[2],label=stat.txt, cex=cex)
    return(mp)
  }
  if(type=="bar"){
    ##browser()
    if(is.null(mar))mar <- c(2.5,3.5,0.5,0.5)
    if(is.null(mgp)) mgp <- c(2.2,0.8,0)
    par(mar=mar, mgp=mgp)
    y.by.x <- table(y,x)
    y.prop <- round(prop.table(y.by.x, margin=2)[2,]*100)
    y.txt <- paste(y.by.x[2,], "/", apply(y.by.x,2, sum), sep="")
    y.lim <- range(0, 105)
    if(is.null(xnames)) xnames <- levels(x)
    if(length(y.prop)==2)
      mp <- barplot(y.prop, space=c(0.3,rep(0.3, length(xnames)-1))/0.6,
                    width=0.6,
                    xlim=c(0,length(xnames)),col="black", names="",
                    ylim=y.lim, xlab="", ylab=y.plab, las=1,...)
    else
      mp <- barplot(y.prop, space=c(0.2,rep(0.4, length(xnames)-1))/0.6,
                    width=0.6,
                    xlim=c(0,length(xnames)),col="black", names="",
                    ylim=y.lim, xlab="", ylab=y.plab, las=1,...)
    axis(1, at=mp, labels=xnames,tick=0,...)
    text(x=mp, y=y.prop+3, y.txt, cex=ifelse(is.null(cex.bar),0.7, cex.bar), font=2)
    if(is.null(loc.stat))
      text(x=mp[1]-0.2*diff(mp)[1], y=98, label=stat.txt, cex=cex, pos=4)
    else
      text(loc.stat[1], loc.stat[2],label=stat.txt, cex=cex)
    axis(2, las=1,...)
    box(...)
    return(mp)
  }
  if(type == "bar2"){
    ##browser()
    if(is.null(mar))mar <- c(2.5,3.5,0.5,0.5)
    if(is.null(mgp)) mgp <- c(2.2,0.8,0)
    par(mar=mar, mgp=mgp)
    y.by.x <- table(y,x)
    y.prop <- round(prop.table(y.by.x, margin = 2)*100)
    y.txt <- paste0(y.by.x, "/", rep(apply(y.by.x, 2, sum),
                                     each = length(levels(y))))
    y.lim <- range(0, 105)
    if(is.null(xnames)) xnames <- levels(x)
    mp <- barplot(y.prop, width=0.6,
                  col= heat.colors(length(levels(y))), names=levels(x),
                  ylim=y.lim, xlab="", ylab=y.plab, las=1, beside = TRUE,
                  legend = levels(y),...)
    text(x=mp, y=y.prop+3, y.txt,
         cex=ifelse(is.null(cex.bar),0.7, cex.bar), font=2)
    if(is.null(loc.stat))
      text(x=mp[1]-0.2*diff(mp)[1], y=98, label=stat.txt, cex=cex, pos=4)
    else
      text(loc.stat[1], loc.stat[2],label=stat.txt, cex=cex)
    axis(2, las=1,...)
    box(...)
    return(mp)
  }
  if(!is.null(fname))  dev.off()
}
