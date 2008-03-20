`siarsinglegrouprun` <-
function(siardata) {
# This function runs the single group MCMC for siar

if(siardata$SHOULDRUN==FALSE || siardata$GRAPHSONLY ==TRUE) {
    cat("You must load in some data first (via option 1) in order to use this feature of the program. \n")
    cat("Press <Enter> to continue")
    readline()
    invisible()
    return(NULL)
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

if(siardata$numgroups == 1) {

# Run size
runchoices <- c("Standard","Long","Very long")
runtitle <- "Choose the size of the model run:"
BADRUN <- TRUE
while(BADRUN ==TRUE) {
    runchoose <- menu(runchoices,title = runtitle)
    if(any(runchoose==seq(1,3))) BADRUN <- FALSE
}

# Now run the code
if(runchoose == 1) output <- siarmcmcdirichlet(siardata$targets,siardata$sources,siardata$corrections)
if(runchoose == 2) output <- siarmcmcdirichlet(siardata$targets,siardata$sources,siardata$corrections,400000,200000,10000,100)
if(runchoose == 3) output <- siarmcmcdirichlet(siardata$targets,siardata$sources,siardata$corrections,1000000,400000,20000,300)

return(output)

} else {

cat("This data has multiple groups - choose the multi group option instead. \n \n")
return(NULL)

}

}

