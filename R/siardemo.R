siardemo <-
function(siarversion=0) {

cat("==================== Demo ==================== \n \n")
cat("This is a simple demo using the Geese data provided \n")
cat("with the package. The data are loaded into the R workspace \n")
cat("but are also available as text files in the package directory. \n")
cat("The data are called geese1demo, sourcesdemo, correctionsdemo \n")
cat("and concdepdemo. \n \n")
cat("This example deals with data with 1 group and 2 isotopes. \n")

geese1demo <- matrix(c(10.22,10.37,10.44,10.52,10.19,10.45,9.91,11.27,9.34,-11.36,-11.88,-10.6,-11.25,-11.66,-10.41,-10.88,-14.73,-11.52),ncol=2,nrow=9)
colnames(geese1demo) <- c("d15NPl","d13CPl")
sourcesdemo <- matrix(c(6.48898447009942,4.4321601032975,11.1926127956759,9.8162797508688,1.45946324319852,2.2680708954525,1.112438464109,0.827103932159427,-11.1702276840033,-30.8798439532492,-11.1709000370306,-14.0570109233264,1.21495615708948,0.641318210364447,1.95933055832669,1.17246767177711),ncol=4,nrow=4)
sourcesdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),sourcesdemo)
colnames(sourcesdemo) <- c("Sources","Meand15N","SDd15N","Meand13C","SDd13C")
correctionsdemo <- matrix(c(rep(3.54,4),rep(0.74,4),rep(1.63,4),rep(0.63,4)),ncol=4,nrow=4)
correctionsdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),correctionsdemo)
colnames(correctionsdemo) <- c("Source","Mean15N","Sd15N","Mean13C","Sd13C")
concdepdemo <- matrix(c(0.0297,0.0355,0.0192,0.0139,0.0097,0.0063,0.0053,0.0057,0.3593,0.4026,0.2098,0.1844,0.0561,0.038,0.0327,0.1131),ncol=4,nrow=4)
concdepdemo <- data.frame(c("Zostera","Grass","U.lactuca","Enteromorpha"),concdepdemo)
colnames(concdepdemo) <- c("Sources","Meand15N","SDd15N","Meand13C","SDd13C")
cat("Press <Enter> to see the data \n")
readline()
invisible()

cat(" The target isotope data is called geese1demo and \n")
cat(" has the following format: \n")
print(geese1demo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The source isotope data is called sourcesdemo and looks like this: \n")
print(sourcesdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The fraction correction data is called correctionsdemo: \n")
print(correctionsdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n The concentration depdendence data is called concdepdemo: \n")
print(concdepdemo)
cat("Press <Enter> to continue... \n")
readline()
invisible()
cat("\n")

cat("The data can be plotted and looks like this... \n")
cat("Press <Enter> to continue \n")
readline()
invisible()

xmins <- min(c(sourcesdemo[,2]-3*sourcesdemo[,3],geese1demo[,1]+correctionsdemo[1,2]))
xmaxs <- max(c(sourcesdemo[,2]+3*sourcesdemo[,3],geese1demo[,1]+correctionsdemo[1,2]))
ymins <- min(c(sourcesdemo[,4]-3*sourcesdemo[,5],geese1demo[,2]+correctionsdemo[2,2]))
ymaxs <- max(c(sourcesdemo[,4]+3*sourcesdemo[,5],geese1demo[,2]+correctionsdemo[2,2]))
newgraphwindow()
plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main="Geese demo data",xlab="d15N",ylab="d13C")
for(j in 1:nrow(geese1demo)) {
    points(geese1demo[j,1]-correctionsdemo[1,2],geese1demo[j,2]-correctionsdemo[2,2],col="grey")
    lines(c(geese1demo[j,1]-correctionsdemo[1,2]-2*correctionsdemo[1,3],geese1demo[j,1]-correctionsdemo[1,2]+2*correctionsdemo[1,3]),c(geese1demo[j,2]-correctionsdemo[2,2],geese1demo[j,2]-correctionsdemo[2,2]),col="grey")
    lines(c(geese1demo[j,1]-correctionsdemo[1,2],geese1demo[j,1]-correctionsdemo[1,2]),c(geese1demo[j,2]-correctionsdemo[2,2]-2*correctionsdemo[2,3],geese1demo[j,2]-correctionsdemo[2,2]+2*correctionsdemo[2,3]),col="grey")
}
points(cbind(sourcesdemo[,2],sourcesdemo[,4]),pch=15,col=seq(1,nrow(sourcesdemo)))
for(i in 1:nrow(sourcesdemo)) {
    lines(c(sourcesdemo[i,2]-2*sourcesdemo[i,3],sourcesdemo[i,2]+2*sourcesdemo[i,3]),c(sourcesdemo[i,4],sourcesdemo[i,4]),col=i)
    lines(c(sourcesdemo[i,2],sourcesdemo[i,2]),c(sourcesdemo[i,4]-2*sourcesdemo[i,5],sourcesdemo[i,4]+2*sourcesdemo[i,5]),col=i)
}
if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
legend(9.5,-24,legend=c(as.character(sourcesdemo[,1]),"data"),lty=c(rep(1,nrow(sourcesdemo)),1),pch=c(rep(15,nrow(sourcesdemo)),1),col=c(seq(1,nrow(sourcesdemo)),"grey"),bty="n")

cat("Press <Enter> to continue \n")
readline()
invisible()


cat("The data can be loaded by typing siarmenu() at the command prompt \n")
cat("and following the options to load in data (option 1), then load \n")
cat("in R objects and follow the instructions. \n \n")
cat("Option 2 runs SIAR for a single group, \n")
cat("This will run for ~10 seconds ... \n")

cat("Press <Enter> to run the model \n")
readline()
invisible()

out <- siarmcmcdirichletv4(geese1demo,sourcesdemo,correctionsdemo,concdepdemo)
out <- out$output

cat("Press <Enter> to continue \n")
readline()
invisible()

cat("From the options menu you can now choose a plot, such as this \n")
cat("density plot... \n")

cat("Press <Enter> to see the plot \n")
readline()
invisible()

sourcenames <- as.character(sourcesdemo[,1])
numsources <- nrow(sourcesdemo)
numdata <- nrow(geese1demo)
numiso <- ncol(sourcesdemo)/2-1
groupnum <- 1
numgroups <- 1
fullname <- "Geese plasma"


usepars <- out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso))]

mybreaks <- seq(0,1,length=50)
halfwidth <- diff(mybreaks)[1]/2
top <- 0
for(j in 1:numsources) {
    top <- max(c(top,max(hist(usepars[,j],plot=FALSE,breaks=mybreaks)$density)))
}
newgraphwindow()
if(fullname!=0) {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities",sep=""),xlab="proportion",ylab="density")
} else {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main="Proportion densities",xlab="proportion",ylab="density")
}
if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:numsources) {
    Ans <- hist(usepars[,j],plot=FALSE,breaks=mybreaks)
    for(k in 1:length(Ans$mids)) {
        lines(c(Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(numsources+1)/2,lend=1)
    }
}
legend(0.6,15,legend=sourcenames,col=seq(1,5),lty=1,lwd=3,bty="n")


cat("Press <Enter> to continue")
readline()
invisible()

cat("With more complicated data sets (see geese2demo), you can fit the model \n")
cat("to multiple groups and produce different types of plots \n")
cat("For advanced users, the function siarmcmcdirichletv4() will allow runs \n")
cat("with different run parameters (such as the number of iterations). \n")
cat("Type help(siarmcmcdirichletv4) for more details. \n \n")

cat("Good luck using the software. \n")
cat("Please report bugs to Andrew.Parnell@ucd.ie \n")

}

