`siarmenu` <-
function() {

library(coda)
library(hdrcde)
library(MASS)

siarversion <-"2.0"

cat("------------------------------- \n")
cat(paste("Welcome to Stable Isotope Analysis in R version", siarversion, "\n"))
cat("Author: Andrew Parnell, Trinity College Dublin\n")
cat("Please report bugs to: Andrew.Parnell@tcd.ie\n")
cat("------------------------------- \n")
cat("\n")
cat("Useful: Press 0 at a prompt to return to the main menu or Esc to exit. \n \n")

siardata <- list()
siardata$SHOULDRUN <- FALSE
siardata$EXIT <- FALSE

while(siardata$EXIT==FALSE)
{

choices <- c("Load in some data (run this first)","Run SIAR for a single group","Run SIAR for multiple groups","Plot single group proportions","Matrix plot of proportions","Plot of proportions by source","Save parameter output to a file","Summary information and convergence diagnostics","Demo (for first time users)","Exit")
title <- "The available options are:"
choose <- menu(choices,title = title)

#####################################################################################################

# Section 1
if(choose == 1) {

siardata <- loadsiardata(siarversion)

cat("Press <Enter> to continue")
readline()
invisible()

}



#####################################################################################################
# Run the single group version

# Section 2
if(choose == 2) {

siardata$output <- siarsinglegrouprun(siardata)

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 3
if(choose == 3) {

siardata$output <- siarmultigrouprun(siardata)

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 4
if(choose == 4) {

siarhistograms(siardata,siarversion)

cat("Press <Enter> to continue")
readline()
invisible()


}
#####################################################################################################

# Section 5
if(choose == 5) {

siarmatrixplot(siardata,siarversion)

cat("Press <Enter> to continue")
readline()
invisible()

}

#####################################################################################################

# Section 6
if(choose == 6) {

siarproportionbysourceplot(siardata,siarversion)

cat("Please maximise this graph before saving or printing. \n")
cat("Press <Enter> to continue")
readline()
invisible()


}

#####################################################################################################

# Section 7
if(choose == 7) {

siarsaveoutput(siardata)

}

#####################################################################################################

if(choose == 8)
{

siarhdrs(siardata)

cat("Press <Enter> to continue...")
readline()
invisible()
    

}

#####################################################################################################

if(choose == 9) {

siardemo(siarversion)

cat("Press <Enter> to continue")
readline()
invisible()

}


#####################################################################################################

if(choose == 10)
{
cat("Thank you. Exiting... \n")
siardata$EXIT=TRUE
}

}

}

