##--------------------------------------------------------
## functions for generating the heatmap
##---------------------------------------------------------
cols2 <- function(low = lowcol, med=medcol, high = highcol, ncolors = 123) {
  ## This function is copied from the heatplot function in library(made4)
  if (is.character(low))
    low <- col2rgb(low)/255
  if (is.character(high))
    high <- col2rgb(high)/255
  if (is.character(med))
    med <- col2rgb(med)/255
  col1 <- rgb(seq(low[1], med[1], len = round(ncolors/2)),
              seq(low[2], med[2], len = round(ncolors/2)),
              seq(low[3], med[3], len = round(ncolors/2)))
  col2 <- rgb(seq(med[1], high[1], len = ncolors-round(ncolors/2)),
              seq(med[2], high[2], len = ncolors-round(ncolors/2)),
              seq(med[3], high[3], len = ncolors-round(ncolors/2)))
  col <- c(col1, col2)
  return(col)
}
fheatmap.gene <- function(dat=dat.work,
                          filename="Figs/fig_heatmap.wmf",
                          lht=c(1,7),
                          width=5, height=7,
                          mar1=c(0,5,1,0.3),
                          mar2=c(3,5,0.2,0.3),
                          lby=1,
                          llen=5,
                          x.axis.line=0,line.top=-0.5,
                          vline.at=NULL, y.cex=0.4){
  ##browser()
  ## heatmaps of biochemicals significantly altered in study samples

  x <- dimnames(dat)[[2]]
  x.ord <- unique(x)
  img.xlab <- gsub("\\.", "\n", x.ord)
  ylab <- dimnames(dat)[[1]]

  ## plotting parameters
  xat.tk <- cumsum(table(x)[x.ord])+0.5
  xat.lab <- cumsum(table(x)[x.ord])-table(x)[x.ord]/2
  xlab <- img.xlab
  Z <- t(dat)

  ## Generate the plot
  if(!is.null(filename))
    win.metafile(file=filename, width=width, height=height)
  layout(matrix(c(1,1, 2,2),nrow=2, ncol=2, byrow=T), height=lht)
  par(mar=mar1)
  image(cbind(seq(1,dim(dat)[2], by=lby)),
        col=cols2(low="green", med="grey10",high="red", ncolors=120), axes=F)
  axis(3, at=seq(0,1, length=llen),
       label=round(seq(1,dim(dat)[2], length=llen)),
       tick=0, line=line.top,cex.axis=0.6)
  par(mar=mar2)
  ## plot the heatmap
  image(x=1:dim(Z)[1], y=1:dim(Z)[2],apply(Z,2,function(x)rank(x)),
        col=cols2(low="green", med="grey10",high="red", ncolors=120), axes=F,
        xlab="", ylab="")
  axis(1, at=c(0,xat.tk), lab=rep("",length(xat.tk)+1))
  axis(1, at=xat.lab,lab=xlab,tick=F, line=x.axis.line, cex.axis=0.8)
  if(is.null(vline.at))vline.at <- seq(2,length(xat.tk), by=2)
  segments(x0=xat.tk[vline.at],
           x1=xat.tk[vline.at],
           y0=rep(0.5, length(vline.at)),
           y1=rep(dim(Z)[2]+0.5, length(vline.at)),
           lwd=1.5, col="white")
  par(mgp=c(2,0.6,0))
  axis(2,at=1:length(ylab),
       label=ylab,
       cex.axis=y.cex, las=1, padj=0)
  if(!is.null(filename)) dev.off()
}
fheatmap.mets <- function(dat=dat.work,
                          out.smry=out.wayne.sgn140619$out.smry,
                          filename="Figs/fig_heatmap2_Wayne_20140620.wmf"){
  ##browser()
  ## heatmaps of biochemicals significantly altered in study samples

  x <- factor(paste(dat[, "Menopausal"],dat[, "BMI.cat2"],dat[,"CLSB"], sep="."),
              levels=c("Pre.Normal.No", "Post.Normal.No",
                       "Pre.Normal.Yes","Post.Normal.Yes",
                       "Pre.OwOb.No","Post.OwOb.No",
                       "Pre.OwOb.Yes","Post.OwOb.Yes"))
  img.xlab <- c("Pre\nNorm\nCLSB-","Post\nNorm\nCLSB-",
                "Pre\nNorm\nCLSB+","Post\nNorm\nCLSB+",
                "Pre\nOwOb\nCLSB-","Post\nOwOb\nCLSB-",
                "Pre\nOwOb\nCLSB+","Post\nOwOb\nCLSB+")
  cid <- order(x)
  ##1. identify the important metabolites by
  ##rid.mets <- which(out.smry[,"OwOb.N.pw.sig"]==-1&
  ##                  out.smry[,"OwOb.N.pt.sig"]==-1&
  ##                  out.smry[,"sig10.bma.CLSB.Int"]!=0)
  rid.mets <- which((out.smry[,"CLSB.pw.sig"]==-1&
                       out.smry[,"CLSB.pt.sig"]==-1)|
                      (out.smry[,"CLSB.pw.sig"]==1&
                         out.smry[,"CLSB.pt.sig"]==1))
  id.img <- rid.mets[order(out.smry[rid.mets,"CLSB.pw"])]
  names.id.img <- as.vector(out.smry[id.img,"WSGNID"])

  ## plotting parameters
  xat.tk <- cumsum(table(x))+0.5
  xat.lab <- cumsum(table(x))-table(x)/2
  xlab <- img.xlab
  xat.lab2 <- NULL
  xlab2 <- NULL
  Z <- dat[cid, names.id.img]

  ## Generate the plot
  if(!is.null(filename))
    win.metafile(file=filename, width=6, height=2.0)
  layout(matrix(c(1,1, 2,2),nrow=2, ncol=2, byrow=T), height=c(1,7))
  par(mar=c(0,5,1,0.3))
  image(cbind(seq(1,dim(dat)[1], length=20)),
        col=cols2(low="green", med="grey10",high="red", ncolors=120), axes=F)
  axis(3, at=seq(0,1, length=4), label=round(seq(1,dim(dat)[1], length=4)),
       tick=0, line=-0.9,cex.axis=0.6)
  par(mar=c(3,5,0.2,0.3))
  ## plot the heatmap
  image(x=1:dim(Z)[1], y=1:dim(Z)[2],apply(Z,2,function(x)rank(x)),
        col=cols2(low="green", med="grey10",high="red", ncolors=120), axes=F,
        xlab="", ylab="")
  axis(1, at=xat.tk, lab=rep("",length(xat.tk)))
  axis(1, at=xat.lab,lab=xlab,tick=F, line=0.5, cex.axis=0.8)
  segments(x0=xat.tk[seq(2,8, by=2)], x1=xat.tk[seq(2,8, by=2)],
           y0=rep(0.5, length(xat.tk)/2),
           y1=rep(dim(Z)[2]+0.5, length(xat.tk)/2), lwd=1.5, col="white")
  ##if(!is.null(xat.lab2))
  ##    axis(1, at=xat.lab2,lab=xlab2,tick=F, line=-0.1, cex=0.8)
  ##if(!is.null(vat)) abline(v=vat, col="gray")
  ##axis(2,at=1:dim(Z)[2], label=out.all[id.img.f,"BIOCHEMICAL"],
  ##     cex.axis=0.6, las=1)
  ylab.col <- as.vector(factor(out.smry[,"LipidType"],
                               lab=c("black", "blue","red")))[id.img]
  unicols <- unique(ylab.col)
  for(i in 1:length(unicols)){
    axis(2,at=grep(unicols[i],ylab.col),
         label=out.smry[id.img,"LipidName"][grep(unicols[i],ylab.col)],
         cex.axis=0.6, las=1,
         col.axis=unicols[i])
  }
  if(!is.null(filename)) dev.off()
}
