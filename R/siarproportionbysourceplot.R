siarproportionbysourceplot <- function(siardata,siarversion=0) {

if(siardata$SHOULDRUN==FALSE && siardata$GRAPHSONLY==FALSE) {
    cat("You must load in some data first (via option 1) in order to use \n")
    cat("this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(NULL)
}

if(length(siardata$output)==0) {
    cat("No output found - check that you have run the SIAR model. \n \n")
    return(NULL)
}

if(siardata$numgroups<2) {
    cat("Number of groups = 1 - cannot run this option \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(NULL)
}

cat("Plot of proportions by source \n")
cat("This requires more than one group in the output file. \n")

cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(length(siardata$sources)>0) {
    sourcenames <- as.character(siardata$sources[,1])
} else {
    sourcenames <- strsplit(colnames(siardata$output[,((groupnum-1)*(siardata$numsources+siardata$numiso)+1):(groupnum*(siardata$numsources+siardata$numiso)-siardata$numiso)]),paste("G",groupnum,sep=""))
}

cat("Enter the source number you wish to plot \n")
cat("The choices are:\n")
title <- "The available options are:"
choose2 <- menu(sourcenames)

# Get the relevant data
groupseq <- seq(1,siardata$numgroups,by=1)
usepars <- siardata$output[,seq(choose2,ncol(siardata$output),by=siardata$numsources+siardata$numiso)]

newgraphwindow()
if(siardata$TITLE!="SIAR data") {
    plot(1,1,xlab="Group",ylab="Proportion",main=paste(siardata$TITLE," by source: ",sourcenames[choose2],sep=""),xlim=c(min(groupseq),max(groupseq)),ylim=c(0,1),type="n",xaxp=c(min(groupseq)-1,max(groupseq)+1,max(groupseq)+1-(min(groupseq)-1)))
} else {
    plot(1,1,xlab="Group",ylab="Proportion",main=paste("Proportions by source: ",sourcenames[choose2],sep=""),xlim=c(min(groupseq),max(groupseq)),ylim=c(0,1),type="n",xaxp=c(min(groupseq)-1,max(groupseq)+1,max(groupseq)+1-(min(groupseq)-1)))
}
if(siarversion>0) mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:ncol(usepars)) {
    temp <- hdr(usepars[,j],c(95,75,50),h=bw.nrd0(usepars[,j]))$hdr
    temp2 <- temp[1,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),lwd=2,lend=2)
    temp2 <- temp[2,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),lwd=6,lend=2)
    temp2 <- temp[3,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),lwd=10,lend=2)
}

legnames <- c("95% error","75% error","50% error")
legend(mean(c(min(groupseq),max(groupseq))),1.02,legend=legnames,lwd=c(2,6,10),ncol=3,xjust=0.5,text.width=strwidth(legnames)/2,bty="n")



}
