siarplotdata <- function(siardata,siarversion=0) {

if(siardata$SHOULDRUN==FALSE) {
  cat("Please enter data via the siarmenu() system of via loadsiardata() \n")
}

# For 2 isotope models
if((ncol(siardata$sources)-1)/2 == 2) {
  cat("\n")
  cat("\n")
  cat("Now plotting data points ...  \n")
  cat("Click to position the legend. \n")
  if(siardata$numgroups==1) {
      xmins <- min(c(siardata$sources[,2]-3*siardata$sources[,3],siardata$targets[,1]-siardata$corrections[1,2]-3*siardata$corrections[1,3]))
      xmaxs <- max(c(siardata$sources[,2]+3*siardata$sources[,3],siardata$targets[,1]-siardata$corrections[1,2]+3*siardata$corrections[1,3]))
      ymins <- min(c(siardata$sources[,4]-3*siardata$sources[,5],siardata$targets[,2]-siardata$corrections[2,2]-3*siardata$corrections[2,3]))
      ymaxs <- max(c(siardata$sources[,4]+3*siardata$sources[,5],siardata$targets[,2]-siardata$corrections[2,2]+3*siardata$corrections[2,3]))
      newgraphwindow()
      if(siardata$corrections[1,1]==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=siardata$TITLE,xlab=colnames(siardata$targets)[1],ylab=colnames(siardata$targets)[2])
      if(siardata$corrections[1,1]!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=siardata$TITLE,xlab=siardata$corrections[1,1],ylab=siardata$corrections[2,1])
      for(j in 1:nrow(siardata$targets)) {
          points(siardata$targets[j,1]-siardata$corrections[1,2],siardata$targets[j,2]-siardata$corrections[2,2],col="grey")
          lines(c(siardata$targets[j,1]-siardata$corrections[1,2]-2*siardata$corrections[1,3],siardata$targets[j,1]-siardata$corrections[1,2]+2*siardata$corrections[1,3]),c(siardata$targets[j,2]-siardata$corrections[2,2],siardata$targets[j,2]-siardata$corrections[2,2]),col="grey")
          lines(c(siardata$targets[j,1]-siardata$corrections[1,2],siardata$targets[j,1]-siardata$corrections[1,2]),c(siardata$targets[j,2]-siardata$corrections[2,2]-2*siardata$corrections[2,3],siardata$targets[j,2]-siardata$corrections[2,2]+2*siardata$corrections[2,3]),col="grey")
      }
      points(cbind(siardata$sources[,2],siardata$sources[,4]),pch=15,col=seq(1,nrow(siardata$sources)))
      for(i in 1:nrow(siardata$sources)) {
          lines(c(siardata$sources[i,2]-2*siardata$sources[i,3],siardata$sources[i,2]+2*siardata$sources[i,3]),c(siardata$sources[i,4],siardata$sources[i,4]),col=i)
          lines(c(siardata$sources[i,2],siardata$sources[i,2]),c(siardata$sources[i,4]-2*siardata$sources[i,5],siardata$sources[i,4]+2*siardata$sources[i,5]),col=i)
      }
      if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
      legend(locator(1),legend=c(as.character(siardata$sources[,1]),"data"),lty=c(rep(1,nrow(siardata$sources)),1),pch=c(rep(15,nrow(siardata$sources)),1),col=c(seq(1,nrow(siardata$sources)),"grey"),bty="n")
  } else {
      xmins <- min(c(siardata$sources[,2]-3*siardata$sources[,3],siardata$targets[,2]-siardata$corrections[1,2]-3*siardata$corrections[1,3]))
      xmaxs <- max(c(siardata$sources[,2]+3*siardata$sources[,3],siardata$targets[,2]-siardata$corrections[1,2]+3*siardata$corrections[1,3]))
      ymins <- min(c(siardata$sources[,4]-3*siardata$sources[,5],siardata$targets[,3]-siardata$corrections[2,2]-3*siardata$corrections[2,3]))
      ymaxs <- max(c(siardata$sources[,4]+3*siardata$sources[,5],siardata$targets[,3]-siardata$corrections[2,2]+3*siardata$corrections[2,3]))
      try(x11(),silent==TRUE)
      if(siardata$corrections[1,1]==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=siardata$TITLE,xlab=colnames(siardata$targets)[2],ylab=colnames(siardata$targets)[3])
      if(siardata$corrections[1,1]!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=siardata$TITLE,xlab=siardata$corrections[1,1],ylab=siardata$corrections[2,1])
      for(j in 1:nrow(siardata$targets)) {
          points(siardata$targets[j,2]-siardata$corrections[1,2],siardata$targets[j,3]-siardata$corrections[2,2],col="grey")
          lines(c(siardata$targets[j,2]-siardata$corrections[1,2]-2*siardata$corrections[1,3],siardata$targets[j,2]-siardata$corrections[1,2]+2*siardata$corrections[1,3]),c(siardata$targets[j,3]-siardata$corrections[2,2],siardata$targets[j,3]-siardata$corrections[2,2]),col="grey")
          lines(c(siardata$targets[j,2]-siardata$corrections[1,2],siardata$targets[j,2]-siardata$corrections[1,2]),c(siardata$targets[j,3]-siardata$corrections[2,2]-2*siardata$corrections[2,3],siardata$targets[j,3]-siardata$corrections[2,2]+2*siardata$corrections[2,3]),col="grey")
      }
      points(cbind(siardata$sources[,2],siardata$sources[,4]),pch=15,col=seq(1,nrow(siardata$sources)))
      for(i in 1:nrow(siardata$sources)) {
          lines(c(siardata$sources[i,2]-2*siardata$sources[i,3],siardata$sources[i,2]+2*siardata$sources[i,3]),c(siardata$sources[i,4],siardata$sources[i,4]),col=i)
          lines(c(siardata$sources[i,2],siardata$sources[i,2]),c(siardata$sources[i,4]-2*siardata$sources[i,5],siardata$sources[i,4]+2*siardata$sources[i,5]),col=i)
      }
      if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
      legend(locator(1),legend=c(as.character(siardata$sources[,1]),"data"),lty=c(rep(1,nrow(siardata$sources)),1),pch=c(rep(15,nrow(siardata$sources)),1),col=c(seq(1,nrow(siardata$sources)),"grey"),bty="n")
  }
}

if((ncol(siardata$sources)-1)/2 > 2) {
cat("Plotting of more than 2 isotopes not supported (yet) \n")
cat("\n")
}

}
