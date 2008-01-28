siarmenu <- function() {

library(coda)
library(hdrcde)

siarversion <-"1.3"

cat("------------------------------- \n")
cat(paste("Welcome to Stable Isotope Analysis in R version", siarversion, "\n"))
cat("Author: Andrew Parnell, Trinity College Dublin\n")
cat("Please report bugs to: Andrew.Parnell@tcd.ie\n")
cat("------------------------------- \n")
cat("\n")
cat("Useful: Press 0 at a prompt to return to the main menu or Esc to exit. \n \n")

PATH <- NULL

SHOULDRUN <- FALSE
GRAPHSONLY <- FALSE
EXIT <- FALSE
while(EXIT==FALSE)
{

choices <- c("Load in some data (run this first)","Run SIAR for a single group","Run SIAR for multiple groups","Plot single group proportions","Matrix plot of proportions","Plot of proportions by source","Save parameter output to a file","Summary information and convergence diagnostics","Demo (for first time users)","Exit")
title <- "The available options are:"
choose <- menu(choices,title = title)

#####################################################################################################

# # Section 1
if(choose == 1) {

choices2 <- c("Load data in from files (Mac not yet supported)","Load in R objects","Load in previous output")
choose2 <- menu(choices2,title = title)

if(choose2==0) siarmenu()

if(choose2==1) {

cat("You need to have created at least 2 text files. \n")
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
    if(PATH==0) siarmenu()
    if(file.exists(PATH)) {
        BADPATH <- FALSE
    } else {
        cat("Cannot find this directory, check your typing \n")
    }
}

BADDATA <- TRUE
while(BADDATA == TRUE) {
    cat("Now input the name of the target isotope file \n")
    cat("(including the file extension eg .txt, .dat, etc) \n")
    DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(DATAFILE)==0) DATAFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(DATAFILE==0) siarmenu()    
    if(file.exists(paste(PATH,"\\",DATAFILE,sep=""))) {
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
    if(SOURCEFILE==0) siarmenu()
    if(file.exists(paste(PATH,"\\",SOURCEFILE,sep=""))) {
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
    if(CORRECTIONSFILE==0) siarmenu()
    if(file.exists(paste(PATH,"\\",CORRECTIONSFILE,sep=""))) {
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
siardata <- as.data.frame(read.table(paste(PATH,"\\",DATAFILE,sep=""),header=TRUE))
sources <- as.data.frame(read.table(paste(PATH,"\\",SOURCEFILE,sep=""),header=TRUE))
if(CORRECTIONSFILE != -999) corrections <- as.data.frame(read.table(paste(PATH,"\\",CORRECTIONSFILE,sep=""),header=TRUE))
cat("Done \n \n")

}

if(choose2==2) {

cat("Please enter the name of the object which contains the target data. \n")
cat("\n")
dataexists <- FALSE
while(dataexists == FALSE) {
    datatemp <- scan(what="",nlines=1,quiet=TRUE)
    while(length(datatemp)==0) datatemp <- scan(what="",nlines=1,quiet=TRUE)
    if(datatemp==0) siarmenu()
    if(!exists(datatemp)) {
        cat("Object not found. Try again or Esc to quit. \n")
    } else {
        siardata <- get(datatemp)
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
    if(sourcestemp==0) siarmenu()
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
        if(correctionstemp==0) siarmenu()
        if(!exists(correctionstemp)) {
            cat("Object not found. Try again or Esc to quit. \n")
        } else {
            corrections <- get(correctionstemp)
            correctionsexists <- TRUE        
        }
    }
}

}

if(choose2==3) {

cat("This option allows you to load in the parameters of a previously saved run \n")
cat("so that you can produce graphs, etc. \n \n")
cat("The structure of this file is such that each row is an iteration and each \n")
cat("column represents either a source proportion parameter or a standard deviation \n")
cat("estimate. More complex analysis can be run by simply loading this file \n")
cat("into R itself (rather than through this program).\n \n")


BADOUTPUT <- TRUE
GRAPHSONLY <- TRUE
while(BADOUTPUT == TRUE) {

    cat("Now input the name of the output file including the file extension \n")
    cat("and the directory where it is located (eg","c:\\siar\\data\\output.txt)",". \n")

    OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    while(length(OUTPUTFILE)==0) OUTPUTFILE <- scan(what="",nlines=1,quiet=TRUE)
    if(OUTPUTFILE==0) siarmenu()
    if(file.exists(OUTPUTFILE)) {
        BADOUTPUT <- FALSE
    } else {
        cat("Cannot find this file, check your typing \n")
    }

}

out <- as.matrix(read.table(file=OUTPUTFILE,header=TRUE))

BADGROUPS <- TRUE
while(BADGROUPS==TRUE) {

    cat("\n")
    cat("Enter the number of groups that this output file represents: \n")
    numgroups <- as.integer(scan(what="",nlines=1,quiet=TRUE))
    while(length(numgroups)==0) numgroups <- scan(what="",nlines=1,quiet=TRUE)
    if(numgroups==0) siarmenu()
    cat("... and the number of isotopes: \n")
    numiso <- as.integer(scan(what="",nlines=1,quiet=TRUE))
    while(length(numiso)==0) numiso <- scan(what="",nlines=1,quiet=TRUE)
    if(numiso==0) siarmenu()
    cat("... and finally the number of sources: \n")
    numsources <- as.integer(scan(what="",nlines=1,quiet=TRUE))
    while(length(numsources)==0) numsources <- scan(what="",nlines=1,quiet=TRUE)
    if(numsources==0) siarmenu()

    if(ncol(out)==numgroups*(numiso+numsources)) {
        cat("Data successfully loaded from previous run. \n")
        cat("You can now plot output and look at summary statistics \n")
        cat("for this particular run. \n") 
        BADGROUPS <- FALSE
    } else {
        cat("The numbers you provided do not match with this output file. \n")
        cat("The output file has",ncol(out),"columns which should be \n")
        cat("the number of groups times (number of isotopes + number of sources) \n")
        cat("Please check your answers \n")
    }

}

}

cat("\n Please enter a name for the data set to be used in plots (or leave blank for default titles). \n")
fullname <- scan(what="",nlines=1,quiet=TRUE,sep="\t")
if(length(fullname)==0) {
    fullname <- 0
} else {
    if(fullname==0) siarmenu()
}

if(choose2==1 || choose2==2) { 
    if(!is.integer(siardata[1,1])) {
        numgroups <- 1
    } else {
        numgroups <- max(siardata[,1])
    }
}

SHOULDRUN <- TRUE

if(GRAPHSONLY == FALSE & ((ncol(sources)-1)/2) == 2) {
cat("\n")
cat("\n")
cat("Now plotting data points ...  \n")
cat("Click to position the legend. \n")
if(corrections[1,1] == 0) corrections <- matrix(0,nrow=1+(ncol(sources)-1)/2,ncol=3)
if(numgroups==1) {
    xmins <- min(c(sources[,2]-3*sources[,3],siardata[,1]-corrections[1,2]-3*corrections[1,3]))
    xmaxs <- max(c(sources[,2]+3*sources[,3],siardata[,1]-corrections[1,2]+3*corrections[1,3]))
    ymins <- min(c(sources[,4]-3*sources[,5],siardata[,2]-corrections[2,2]-3*corrections[2,3]))
    ymaxs <- max(c(sources[,4]+3*sources[,5],siardata[,2]-corrections[2,2]+3*corrections[2,3]))
    windows()
    if(fullname==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main="SIAR data",xlab=corrections[1,1],ylab=corrections[2,1])
    if(fullname!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=fullname,xlab=corrections[1,1],ylab=corrections[2,1])
    for(j in 1:nrow(siardata)) {
        points(siardata[j,1]-corrections[1,2],siardata[j,2]-corrections[2,2],col="grey")
        lines(c(siardata[j,1]-corrections[1,2]-2*corrections[1,3],siardata[j,1]-corrections[1,2]+2*corrections[1,3]),c(siardata[j,2]-corrections[2,2],siardata[j,2]-corrections[2,2]),col="grey")
        lines(c(siardata[j,1]-corrections[1,2],siardata[j,1]-corrections[1,2]),c(siardata[j,2]-corrections[2,2]-2*corrections[2,3],siardata[j,2]-corrections[2,2]+2*corrections[2,3]),col="grey")
    }
    points(cbind(sources[,2],sources[,4]),pch=15,col=seq(1,nrow(sources)))
    for(i in 1:nrow(sources)) {
        lines(c(sources[i,2]-2*sources[i,3],sources[i,2]+2*sources[i,3]),c(sources[i,4],sources[i,4]),col=i)
        lines(c(sources[i,2],sources[i,2]),c(sources[i,4]-2*sources[i,5],sources[i,4]+2*sources[i,5]),col=i)
    }
    mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
    legend(locator(1),legend=c(as.character(sources[,1]),"data"),lty=c(rep(1,nrow(sources)),1),pch=c(rep(15,nrow(sources)),1),col=c(seq(1,nrow(sources)),"grey"),bty="n")
} else {
    xmins <- min(c(sources[,2]-3*sources[,3],siardata[,2]-corrections[1,2]-3*corrections[1,3]))
    xmaxs <- max(c(sources[,2]+3*sources[,3],siardata[,2]-corrections[1,2]+3*corrections[1,3]))
    ymins <- min(c(sources[,4]-3*sources[,5],siardata[,3]-corrections[2,2]-3*corrections[2,3]))
    ymaxs <- max(c(sources[,4]+3*sources[,5],siardata[,3]-corrections[2,2]+3*corrections[2,3]))
    windows()
    if(fullname==0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main="SIAR data",xlab=corrections[1,1],ylab=corrections[2,1])
    if(fullname!=0) plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main=fullname,xlab=corrections[1,1],ylab=corrections[2,1])
    for(j in 1:nrow(siardata)) {
        points(siardata[j,2]-corrections[1,2],siardata[j,3]-corrections[2,2],col="grey")
        lines(c(siardata[j,2]-corrections[1,2]-2*corrections[1,3],siardata[j,2]-corrections[1,2]+2*corrections[1,3]),c(siardata[j,3]-corrections[2,2],siardata[j,3]-corrections[2,2]),col="grey")
        lines(c(siardata[j,2]-corrections[1,2],siardata[j,2]-corrections[1,2]),c(siardata[j,3]-corrections[2,2]-2*corrections[2,3],siardata[j,3]-corrections[2,2]+2*corrections[2,3]),col="grey")
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

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################



# Section 2
if(choose == 2) {

if(SHOULDRUN==FALSE || GRAPHSONLY ==TRUE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Run SIAR for a single group. \n")
cat("For this you will need to have input successfully some data in option 1. \n")
cat("In this instance, the target isotope file must be a two column file with \n")
cat("no group label in the first column. \n \n")
cat("See the demo for more details on this function. \n")
cat("\n")
cat("Press <Enter> to continue...")
readline()
invisible()

# Run size
runchoices <- c("Standard","Long","Very long")
runtitle <- "Choose the size of the model run:"
runchoose <- menu(runchoices,title = runtitle)

# Now run the code
if(runchoose == 1) out <- siarmcmc(siardata,sources,corrections)
if(runchoose == 2) out <- siarmcmc(siardata,sources,corrections,400000,200000,10000,100)
if(runchoose == 3) out <- siarmcmc(siardata,sources,corrections,1000000,400000,20000,300)

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 3
if(choose == 3) {

if(SHOULDRUN==FALSE || GRAPHSONLY==TRUE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Run SIAR for multiple groups. \n")
cat("For this you will need to have input successfully some data in option 1. \n")
cat("In this instance, the target isotope file must be a three column file with \n")
cat("a group label in the first column. \n \n")
cat("See the demo for more details on this function. \n")
cat("\n")
cat("Press <Enter> to continue...")
readline()
invisible()

# Run size
runchoices <- c("Standard","Long","Very long")
runtitle <- "Choose the size of the model run:"
runchoose <- menu(runchoices,title = runtitle)

# Now run the code
if(runchoose == 1) out <- siarmcmc(siardata,sources,corrections)
if(runchoose == 2) out <- siarmcmc(siardata,sources,corrections,400000,200000,10000,20)
if(runchoose == 3) out <- siarmcmc(siardata,sources,corrections,1000000,400000,20000,60)

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 4
if(choose == 4) {

if(SHOULDRUN==FALSE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Plots of single groups proportions. \n")

cat("Enter the group number of the proportions you wish to plot \n")
cat("or leave blank if you have only have 1 group \n")
cat("Click on the graph to position the legend... \n")
groupnum <- as.integer(scan(what="integer",nlines=1,quiet=TRUE))

if(length(groupnum)==0) {
    groupnum <- 1
} else {
    if(groupnum==0) siarmenu()
}

if(!exists("out")) cat("Cannot find output file - run option 1 again \n")

cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(exists("sources")) {
    sourcenames <- as.character(sources[,1])
    numsources <- nrow(sources)
    numiso <- (ncol(sources)-1)/2
} else {
    #sourcenames <- colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)])
    sourcenames <- strsplit(colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)]),paste("G",groupnum,sep=""))
}

# Get the right dimensions of pars
usepars <- out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso))]

mybreaks <- seq(0,1,length=50)
halfwidth <- diff(mybreaks)[1]/2
top <- 0
for(j in 1:numsources) {
    top <- max(c(top,max(hist(usepars[,j],plot=FALSE,breaks=mybreaks)$density)))
}

windows()
if(fullname!=0) {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities for group ",groupnum,sep=""),xlab="roportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities",sep=""),xlab="roportion",ylab="density")
} else {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main="Proportion densities",xlab="proportion",ylab="density")
}
mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:numsources) {
    Ans <- hist(usepars[,j],plot=FALSE,breaks=mybreaks)
    for(k in 1:length(Ans$mids)) {
        lines(c(Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(numsources+1)/2,lend=1)
    }
}
legend(locator(1),legend=sourcenames,col=seq(1,5),lty=1,lwd=3,bty="n")

cat("Press <Enter> to continue")
readline()
invisible()


}
#####################################################################################################

# Section 5
if(choose == 5) {

if(SHOULDRUN==FALSE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Matrix plot of groups proportions. \n")

cat("Enter the group number of the proportions you wish to plot \n")
cat("or leave blank if you have only have 1 group \n")
groupnum <- as.integer(scan(what="integer",nlines=1,quiet=TRUE))
if(length(groupnum)==0) {
    groupnum <- 1
} else {
    if(groupnum==0) siarmenu()
}

if(!exists("out")) {
    cat("Cannot find output file - run options 1 or 2 again \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}


cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(exists("sources")) {
    sourcenames <- as.character(sources[,1])
    numsources <- nrow(sources)
    numiso <- (ncol(sources)-1)/2
} else {
    #sourcenames <- colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)])
    sourcenames <- strsplit(colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)]),paste("G",groupnum,sep=""))
}

# Get some column names
if(exists("siardata")) {
    if(!is.integer(siardata[1,1])) {
        colnames(out) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
    } else {
        colnames(out) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,numiso),sep="")),"G",sep=""),times=numgroups),sort(rep(seq(1,numgroups),times=numsources+numiso)),sep="")     
    }
}

# Get the right dimensions of pars
usepars <- out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso))]

windows()
if(fullname!=0) {
    if(numgroups > 1) pairs(usepars[,1:numsources],xlim=c(0,1),ylim=c(0,1),main=paste(fullname,": matrix plot of proportions for group ",groupnum,sep=""),diag.panel=panelhist,lower.panel=panelcor)
    if(numgroups ==1) pairs(usepars[,1:numsources],xlim=c(0,1),ylim=c(0,1),main=paste(fullname,": matrix plot of proportions",sep=""),diag.panel=panelhist,lower.panel=panelcor)
} else {
    if(numgroups > 1) pairs(usepars[,1:numsources],xlim=c(0,1),ylim=c(0,1),main=paste("Matrix plot of proportions for group ",groupnum,sep=""),diag.panel=panelhist,lower.panel=panelcor)
    if(numgroups ==1) pairs(usepars[,1:numsources],xlim=c(0,1),ylim=c(0,1),main="Matrix plot of proportions",diag.panel=panelhist,lower.panel=panelcor)
}

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 6
if(choose == 6) {

if(SHOULDRUN==FALSE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Plot of proportions by group. \n")
cat("This requires more than one group in the output file. \n")

if(!exists("out")) {
    cat("Cannot find output file - run option 1 again \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}
if(numgroups<2) {
    cat("Number of groups = 1 - cannot run this option \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Producing plot..... \n \n")

# Define some of the useful things the function needs to know
if(exists("sources")) {
    sourcenames <- as.character(sources[,1])
    numsources <- nrow(sources)
    numiso <- (ncol(sources)-1)/2
} else {
    #sourcenames <- colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)])
    sourcenames <- strsplit(colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)]),paste("G",groupnum,sep=""))
}


cat("Enter the source number you wish to plot \n")
cat("The choices are:\n")
title <- "The available options are:"
choose2 <- menu(sourcenames)

# Define some of the useful things the function needs to know
if(exists("sources")) {
    sourcenames <- as.character(sources[,1])
    numsources <- nrow(sources)
    numiso <- (ncol(sources)-1)/2
} else {
    sourcenames <- colnames(out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso)-numiso)])
}

# Get some column names
if(exists("siardata")) {
    if(!is.integer(siardata[1,1])) {
        colnames(out) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
    } else {
        colnames(out) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,numiso),sep="")),"G",sep=""),times=numgroups),sort(rep(seq(1,numgroups),times=numsources+numiso)),sep="")
    }
}

# Get the relevant data
groupseq <- seq(1,numgroups,by=1)
usepars <- out[,seq(choose2,ncol(out),by=numsources+numiso)]

windows()
if(fullname!=0) {
    plot(1,1,xlab="Group",ylab="Proportion",main=paste(fullname," by group: ",sourcenames[choose2],sep=""),xlim=c(min(groupseq),max(groupseq)),ylim=c(0,1),type="n",xaxp=c(min(groupseq)-1,max(groupseq)+1,max(groupseq)+1-(min(groupseq)-1)))
} else {
    plot(1,1,xlab="Group",ylab="Proportion",main=paste("Proportions by group: ",sourcenames[choose2],sep=""),xlim=c(min(groupseq),max(groupseq)),ylim=c(0,1),type="n",xaxp=c(min(groupseq)-1,max(groupseq)+1,max(groupseq)+1-(min(groupseq)-1)))
}
mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:ncol(usepars)) {
    temp <- hdr(usepars[,j],c(95,75,50),h=bw.nrd0(usepars[,j]))$hdr
    temp2 <- temp[1,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),col="yellow",lwd=3)
    temp2 <- temp[2,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),col="orange",lwd=3)
    temp2 <- temp[3,]
    lines(c(groupseq[j],groupseq[j]),c(min(temp2[!is.na(temp2)]),max(temp2[!is.na(temp2)])),col="red",lwd=3)
}

legnames <- c("95% error","75% error","50% error")
legend(mean(c(min(groupseq),max(groupseq))),1.02,legend=legnames,col=c("yellow","orange","red"),lwd=c(3,3,3),ncol=3,xjust=0.5,text.width=strwidth(legnames)/2,bty="n")

cat("Please maximise this graph before saving or printing. \n")
cat("Press <Enter> to continue")
readline()
invisible()


}

#####################################################################################################

# Section 7
if(choose == 7) {

if(SHOULDRUN==FALSE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

BADFILE <- TRUE
while(BADFILE == TRUE) {
    cat("To save the output you need to have run either option 2 or 3 \n")
    cat("from the main menu. \n \n")

    cat("Enter a directory location where the output parameters will reside: \n")

    outputfileloc <- scan(what="",nlines=1,quiet=TRUE)
    while(length(outputfileloc)==0) outputfileloc <- scan(what="",nlines=1,quiet=TRUE)
    if(outputfileloc==0) siarmenu()

    if(!file.exists(outputfileloc)) {
        cat("This location doesn't exist, check your typing \n")
    } else {
        cat("Please enter a filename and a file extension for the parameters: \n")
        outputfilename <- scan(what="",nlines=1,quiet=TRUE)
        while(length(outputfilename)==0) outputfilename <- scan(what="",nlines=1,quiet=TRUE)
        if(outputfilename==0) siarmenu()

        # Define some of the useful things the function needs to know
        if(!exists("siardata")) {
            cat("Error: no data entered, go back and run the model first \n")
            cat("Press <Enter> to continue")
        readline()
        invisible()
        siarmenu()
        }

        sourcenames <- as.character(sources[,1])
        numsources <- nrow(sources)
        numdata <- nrow(siardata)
        numiso <- (ncol(sources)-1)/2

        # Get some column names
        if(!is.integer(siardata[1,1])) {
            colnames(out) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
        } else {
            colnames(out) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,numiso),sep="")),"G",sep=""),times=numgroups),sort(rep(seq(1,numgroups),times=numsources+numiso)),sep="")     
        }


        cat("Writing output ... \n")
        write.table(out,file=paste(outputfileloc,"\\",outputfilename,sep=""),quote=FALSE,row.names=FALSE)
        BADFILE <- FALSE
        cat("Output created. \n \n ")
        
        cat("Press <Enter> to continue")
        readline()
        invisible()

    }

}

}

#####################################################################################################

if(choose == 8)
{

if(SHOULDRUN==FALSE) {
    cat("You must run option 1 or option 9 first in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}

cat("Summary information for the output file ... \n")

if(!exists("out")) {
    cat("Cannot find any data runs - please run option 1 from the main menu. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    siarmenu()
}
hdrsummary <- matrix(0,ncol=4,nrow=ncol(out))
colnames(hdrsummary) <- c("Low 95% hdr","High 95% hdr","mode","mean")
if(exists("siardata")) {
    sourcenames <- as.character(sources[,1])
    numiso <- (ncol(sources)-1)/2
    if(!is.integer(siardata[1,1])) {    
        rownames(hdrsummary) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
    } else {
        rownames(hdrsummary) <- paste(rep(paste(c(sourcenames,paste("SD",seq(1,numiso),sep="")),"G",sep=""),times=numgroups),sort(rep(seq(1,numgroups),times=numsources+numiso)),sep="")
    }
} else {
    rownames(hdrsummary) <- colnames(out)
}

for(i in 1:ncol(out)) {
    temp <- hdr(out[,i],h=bw.nrd0(out[,i]))
    hdrsummary[i,1] <- max(0,temp$hdr[2,1])
    hdrsummary[i,2] <- temp$hdr[2,2]
    hdrsummary[i,3] <- temp$mode
    hdrsummary[i,4] <- mean(out[,i])
}

print(hdrsummary)

cat("Press <Enter> to continue...")
readline()
invisible()

cat("\n")
cat("Running convergence diagnostics on output. \n")
cat("Output parameters need to have been loaded in or created. \n \n")

cat("Worst parameters are ... \n")
temp <- geweke.diag(out)[[1]]
print(sort(c(pnorm(temp[temp<0]),1-pnorm(temp[temp>0])))[1:min(10,ncol(out))])

cat("If lots of the p-values are very small , try a longer run of the MCMC. \n")

cat("Press <Enter> to continue...")
readline()
invisible()
    

}

#####################################################################################################

if(choose == 9) {

cat("==================== Demo ==================== \n \n")
cat("This is a simple demo using the Geese data provided \n")
cat("with the package. The data are loaded into the R workspace \n")
cat("but are also available as text files in the package directory. \n")
cat("The data are called geese1demo, sourcesdemo and correctionsdemo. \n \n")
cat("This example deals with data with 1 group and 2 isotopes. \n")

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
cat("\n")

cat("The data can be plotted and looks like this... \n")
cat("Press <Enter> to continue \n")
readline()
invisible()

xmins <- min(c(sourcesdemo[,2]-3*sourcesdemo[,3],geese1demo[,1]+correctionsdemo[1,2]))
xmaxs <- max(c(sourcesdemo[,2]+3*sourcesdemo[,3],geese1demo[,1]+correctionsdemo[1,2]))
ymins <- min(c(sourcesdemo[,4]-3*sourcesdemo[,5],geese1demo[,2]+correctionsdemo[2,2]))
ymaxs <- max(c(sourcesdemo[,4]+3*sourcesdemo[,5],geese1demo[,2]+correctionsdemo[2,2]))
windows()
plot(1,1,type="n",xlim=c(xmins,xmaxs),ylim=c(ymins,ymaxs),main="Geese demo data",xlab=correctionsdemo[1,1],ylab=correctionsdemo[2,1])
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
mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)
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

out <- siarmcmc(geese1demo,sourcesdemo,correctionsdemo)

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

# Get some column names
colnames(out) <- c(sourcenames,paste("SD",seq(1,numiso),sep=""))
  
usepars <- out[,((groupnum-1)*(numsources+numiso)+1):(groupnum*(numsources+numiso))]

mybreaks <- seq(0,1,length=50)
halfwidth <- diff(mybreaks)[1]/2
top <- 0
for(j in 1:numsources) {
    top <- max(c(top,max(hist(usepars[,j],plot=FALSE,breaks=mybreaks)$density)))
}
windows()
if(fullname!=0) {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste(fullname,": proportion densities",sep=""),xlab="proportion",ylab="density")
} else {
    if(numgroups > 1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main=paste("Proportion densities for group ",groupnum,sep=""),xlab="proportion",ylab="density")
    if(numgroups ==1) plot(1,1,xlim=c(0,1),ylim=c(0,top),type="n",main="Proportion densities",xlab="proportion",ylab="density")
}
mtext(paste("siar v",siarversion),side=1,line=4,adj=1,cex=0.6)

for(j in 1:numsources) {
    Ans <- hist(usepars[,j],plot=FALSE,breaks=mybreaks)
    for(k in 1:length(Ans$mids)) {
        lines(c(Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth,Ans$mids[k]+(j/((numsources+1)/2)-1)*halfwidth),c(0,Ans$density[k]),col=j,lwd=(numsources+1)/2,lend=1)
    }
}
legend(0.6,12,legend=sourcenames,col=seq(1,5),lty=1,lwd=3,bty="n")


cat("Press <Enter> to continue")
readline()
invisible()

cat("With more complicated data sets (see geese2demo), you can fit the model \n")
cat("to multiple groups and produce different types of plots \n")
cat("For advanced users, the function siarmcmc() will allow runs \n")
cat("with different run parameters (such as the number of iterations). \n")
cat("Type help(siarmcmc) for more details. \n \n")

cat("Good luck using the software. \n")
cat("Please report bugs to Andrew.Parnell@tcd.ie \n")

cat("Press <Enter> to continue")
readline()
invisible()

}


#####################################################################################################

if(choose == 10)
{
cat("Thank you. Exiting... \n")
EXIT=TRUE
}

}

}
