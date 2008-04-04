loadsiardata <- function(siarversion) {

choices2 <- c("Load data in from files","Load in R objects","Load in previous output")
title <- "The available options are:"
choose2 <- menu(choices2,title = title)

############################################################################

if(choose2==0) return(list(EXIT=FALSE))

############################################################################
# Load in data from files

if(choose2==1) {

cat("To run siar, you need to have created at least 2 text files. \n")
cat("The first must contain the target isotope measurements in either \n")
cat("two columns with no group number or 3 columns with a group label \n")
cat("The second file must contain a column of the different source names \n")
cat("followed by the isotope measurements for each in a seperate columns \n \n")
cat("Optionally, a third file can be created which contains the fractionation \n")
cat("correction means and standard deviations for each isotope \n \n")
cat("See the demo and the included data files for more information \n")
cat("On the data input format \n \n")

BADPATH <- TRUE
while(BADPATH == TRUE) {
    cat("First input the directory at which the files can be found: \n")
    PATH <- scan(what="",nlines=1,quiet=TRUE)
    while(length(PATH)==0) PATH <- scan(what="",nlines=1,quiet=TRUE)
    if(PATH==0) return(list(EXIT=FALSE))
    if(file.exists(PATH)) {
        BADPATH <- FALSE
    } else {
        cat("Cannot find this directory. Check your typing. \n")
    }
}

BADDATA <- TRUE
while(BADDATA == TRUE) {
    cat("Now input the name of the target isotope file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(DATAFILE)==0) DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(DATAFILE==0) return(list(EXIT=FALSE))    
    if(file.exists(paste(PATH,"/",DATAFILE,sep=""))) {
        BADDATA <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }
}

BADSOURCES <- TRUE
while(BADSOURCES == TRUE) {
    cat("Now input the name of the source isotope file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")    
    SOURCEFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(SOURCEFILE)==0) SOURCEFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(SOURCEFILE==0) return(list(EXIT=FALSE))
    if(file.exists(paste(PATH,"/",SOURCEFILE,sep=""))) {
        BADSOURCES <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }
}

BADCORRECTIONS <- TRUE
while(BADCORRECTIONS == TRUE) {
    cat("Now input the name of the fractionation corrections file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    cat("or leave blank to use pre-corrected values \n")
    CORRECTIONSFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(length(CORRECTIONSFILE)==0) CORRECTIONSFILE <- -999
    if(CORRECTIONSFILE==0) return(list(EXIT=FALSE))
    if(file.exists(paste(PATH,"/",CORRECTIONSFILE,sep=""))) {
        BADCORRECTIONS <- FALSE
    } else {
        if(CORRECTIONSFILE == -999) {
            BADCORRECTIONS <- FALSE
            corrections <- matrix(0,nrow=1,ncol=1)
        }
        if(CORRECTIONSFILE !=-999) cat("Cannot find this file, check your typing \n")
    }
}

cat("Now loading in data... \n")
targets <- as.data.frame(read.table(paste(PATH,"/",DATAFILE,sep=""),header=TRUE))
sources <- as.data.frame(read.table(paste(PATH,"/",SOURCEFILE,sep=""),header=TRUE))
if(CORRECTIONSFILE != -999) corrections <- as.data.frame(read.table(paste(PATH,"/",CORRECTIONSFILE,sep=""),header=TRUE))
cat("Done \n \n")

# Finally sort everything out so its in proper siar format
if(corrections[1,1] == 0) corrections <- matrix(0,nrow=1+(ncol(sources)-1)/2,ncol=3)
numgroups <- 1
if(targets[1,1]%%1 == 0) numgroups <- max(targets[,1])
numsources <-nrow(sources)
numdata <- nrow(targets)
numiso <- (ncol(sources)-1)/2
SHOULDRUN <- TRUE
GRAPHSONLY <- FALSE
EXIT <- FALSE
output <- NULL

}

############################################################################
# Load in R objects

if(choose2==2) {

cat("Please enter the name of the object which contains the target data. \n")
cat("\n")
dataexists <- FALSE
while(dataexists == FALSE) {
    datatemp <- scan(what="",nlines=1,quiet=TRUE)
    while(length(datatemp)==0) datatemp <- scan(what="",nlines=1,quiet=TRUE)
    if(datatemp==0) return(list(EXIT=FALSE))
    if(!exists(datatemp)) {
        cat("Object not found. Try again or Esc to quit. \n")
    } else {
        targets <- get(datatemp)
        dataexists <- TRUE        
    }
}

cat("Now please enter the name of the object which contains the source \n")
cat("isotope details. The first column should be the source names \n")
cat("\n")
sourcesexists <- FALSE
while(sourcesexists == FALSE) {
    sourcestemp <- scan(what="",nlines=1,quiet=TRUE)
    while(length(sourcestemp)==0) sourcestemp <- scan(what="",nlines=1,quiet=TRUE)
    if(sourcestemp==0) return(list(EXIT=FALSE))
    if(!exists(sourcestemp)) {
        cat("Object not found. Try again or Esc to quit. \n")
    } else {
        sources <- get(sourcestemp)
        sourcesexists <- TRUE        
    }
}

cat("Now please enter the name of the object which contains the isotopic \n")
cat("correction mean and standard deviation. Note: if the data are \n")
cat("pre-corrected please leave blank \n")
cat("\n")
correctionsexists <- FALSE
while(correctionsexists == FALSE) {
    correctionstemp <- scan(what="",nlines=1,quiet=TRUE)
    if(length(correctionstemp)==0) {
        corrections <- matrix(0,nrow=1,ncol=1)
        correctionsexists <- TRUE
    } else {
        if(correctionstemp==0) return(list(EXIT=FALSE))
        if(!exists(correctionstemp)) {
            cat("Object not found. Try again or Esc to quit. \n")
        } else {
            corrections <- get(correctionstemp)
            correctionsexists <- TRUE        
        }
    }
}

# Finally sort everything out so its in proper siar format
if(corrections[1,1] == 0) corrections <- matrix(0,nrow=1+(ncol(sources)-1)/2,ncol=3)
numgroups <- 1
if(targets[1,1]%%1 == 0) numgroups <- max(targets[,1])
numsources <-nrow(sources)
numdata <- nrow(targets)
numiso <- (ncol(sources)-1)/2
PATH <- NULL
SHOULDRUN <- TRUE
GRAPHSONLY <- FALSE
EXIT <- FALSE
output <- NULL

}

############################################################################
# Load in a previous run

if(choose2==3) {

cat("This option allows you to load in the parameters of a previously saved run \n")
cat("so that you can produce graphs, etc. \n \n")
cat("More complex analysis can be run by simply loading this file \n")
cat("into R itself via load(data) (rather than through this menu system).\n \n")

BADOUTPUT <- TRUE
GRAPHSONLY <- TRUE
while(BADOUTPUT == TRUE) {

    cat("Now input the name of the output file including the file extension \n")
    cat("and the directory where it is located (eg","c:\\siar\\data\\output.Rdata)",". \n")

    OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(OUTPUTFILE)==0) OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(OUTPUTFILE==0) return(list(EXIT=FALSE))
    if(file.exists(OUTPUTFILE)) {
        BADOUTPUT <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }

}

#output <- as.matrix(read.table(file=OUTPUTFILE,header=TRUE))
load(file=OUTPUTFILE)

# Finally sort everything out so its in proper siar format
targets <- siardata$targets
sources <- siardata$sources
corrections <- siardata$corrections
PATH <- siardata$PATH
numdata <- siardata$numdata
SHOULDRUN <- siardata$SHOULDRUN
GRAPHSONLY <- siardata$GRAPHSONLY
EXIT <- siardata$EXIT
output <- siardata$output
TITLE <- siardata$TITLE
numgroups <- siardata$numgroups
numdata <- siardata$numdata
numsources <- siardata$numsources
numiso <- siardata$numiso

}

############################################################################
# Enter a title

if(choose2==1 || choose2==2) {
    cat("\n Please enter a name for the data set to be used in plots (or leave blank for default titles). \n")
    TITLE <- scan(what="",nlines=1,quiet=TRUE,sep="\t")
    if(length(TITLE)==0) {
        TITLE <- "SIAR data"
    } else {
        if(TITLE==0) return(NULL)
    }
}

############################################################################
# Graphs for 2 isotope version

if(GRAPHSONLY == FALSE) {

BADPLOTGRAPH <- FALSE
while(BADPLOTGRAPH == FALSE) {
    cat("Do you wish to plot the data? (y/n) \n")
    PLOTGRAPH <- scan(what="",nlines=1,quiet=TRUE)
    if(length(PLOTGRAPH)>0) BADPLOTGRAPH <- TRUE
}

if(PLOTGRAPH == "y" || PLOTGRAPH == "yes") {
if((ncol(sources)-1)/2 == 2) {
cat("\n")
cat("\n")
cat("Now plotting data points ...  \n")
cat("Click to position the legend. \n")
if(numgroups==1) {
    xmins <- min(c(sources[,2]-3*sources[,3],targets[,1]-corrections[1,2]-3*corrections[1,3]))
    xmaxs <- max(c(sources[,2]+3*sources[,3],targets[,1]-corrections[1,2]+3*corrections[1,3]))
    ymins <- min(c(sources[,4]-3*sources[,5],targets[,2]-corrections[2,2]-3*corrections[2,3]))
    ymaxs <- max(c(sources[,4]+3*sources[,5],targets[,2]-corrections[2,2]+3*corrections[2,3]))
    newgraphwindow()
    if(corrections[1,1]==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=TITLE,xlab=colnames(targets)[1],ylab=colnames(targets)[2])
    if(corrections[1,1]!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=TITLE,xlab=corrections[1,1],ylab=corrections[2,1])
    for(j in 1:nrow(targets)) {
        points(targets[j,1]-corrections[1,2],targets[j,2]-corrections[2,2],col="grey")
        lines(c(targets[j,1]-corrections[1,2]-2*corrections[1,3],targets[j,1]-corrections[1,2]+2*corrections[1,3]),c(targets[j,2]-corrections[2,2],targets[j,2]-corrections[2,2]),col="grey")
        lines(c(targets[j,1]-corrections[1,2],targets[j,1]-corrections[1,2]),c(targets[j,2]-corrections[2,2]-2*corrections[2,3],targets[j,2]-corrections[2,2]+2*corrections[2,3]),col="grey")
    }
    points(cbind(sources[,2],sources[,4]),pch=15,col=seq(1,nrow(sources)))
    for(i in 1:nrow(sources)) {
        lines(c(sources[i,2]-2*sources[i,3],sources[i,2]+2*sources[i,3]),c(sources[i,4],sources[i,4]),col=i)
        lines(c(sources[i,2],sources[i,2]),c(sources[i,4]-2*sources[i,5],sources[i,4]+2*sources[i,5]),col=i)
    }
    mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
    legend(locator(1),legend=c(as.character(sources[,1]),"data"),lty=c(rep(1,nrow(sources)),1),pch=c(rep(15,nrow(sources)),1),col=c(seq(1,nrow(sources)),"grey"),bty="n")
} else {
    xmins <- min(c(sources[,2]-3*sources[,3],targets[,2]-corrections[1,2]-3*corrections[1,3]))
    xmaxs <- max(c(sources[,2]+3*sources[,3],targets[,2]-corrections[1,2]+3*corrections[1,3]))
    ymins <- min(c(sources[,4]-3*sources[,5],targets[,3]-corrections[2,2]-3*corrections[2,3]))
    ymaxs <- max(c(sources[,4]+3*sources[,5],targets[,3]-corrections[2,2]+3*corrections[2,3]))
    try(x11(),silent==TRUE)
    if(corrections[1,1]==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=TITLE,xlab=colnames(targets)[2],ylab=colnames(targets)[3])
    if(corrections[1,1]!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=TITLE,xlab=corrections[1,1],ylab=corrections[2,1])
    for(j in 1:nrow(targets)) {
        points(targets[j,2]-corrections[1,2],targets[j,3]-corrections[2,2],col="grey")
        lines(c(targets[j,2]-corrections[1,2]-2*corrections[1,3],targets[j,2]-corrections[1,2]+2*corrections[1,3]),c(targets[j,3]-corrections[2,2],targets[j,3]-corrections[2,2]),col="grey")
        lines(c(targets[j,2]-corrections[1,2],targets[j,2]-corrections[1,2]),c(targets[j,3]-corrections[2,2]-2*corrections[2,3],targets[j,3]-corrections[2,2]+2*corrections[2,3]),col="grey")
    }
    points(cbind(sources[,2],sources[,4]),pch=15,col=seq(1,nrow(sources)))
    for(i in 1:nrow(sources)) {
        lines(c(sources[i,2]-2*sources[i,3],sources[i,2]+2*sources[i,3]),c(sources[i,4],sources[i,4]),col=i)
        lines(c(sources[i,2],sources[i,2]),c(sources[i,4]-2*sources[i,5],sources[i,4]+2*sources[i,5]),col=i)
    }
    mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
    legend(locator(1),legend=c(as.character(sources[,1]),"data"),lty=c(rep(1,nrow(sources)),1),pch=c(rep(15,nrow(sources)),1),col=c(seq(1,nrow(sources)),"grey"),bty="n")
}
}
}
}

return(list(targets=targets,sources=sources,corrections=corrections,PATH=PATH,TITLE=TITLE,numgroups=numgroups,numdata=numdata,numsources=numsources,numiso=numiso,SHOULDRUN=SHOULDRUN,GRAPHSONLY=GRAPHSONLY,EXIT=EXIT,output=output))

}
