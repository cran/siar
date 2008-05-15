siarplotdata <- function(siardata, siarversion = 0, grp=1:siardata$numgroups,panel=NULL,isos=c(0,0),leg=1){

    #  SUB FUNCTION TO PLOT THE CROSSES ON THE GRAPH
    siaraddcross <- function (x=NULL,ex=NULL,y=NULL,ey=NULL,clr="grey50",upch=21) {

    points(x, y, col=clr, pch = upch)
    
    if(!is.null(ex)){ lines(c(x-ex, x+ex),c(y,y), col=clr)}
    
    if (!is.null(ey)){ lines(c(x,x), c(y-ey,y+ey), col=clr)}



    }
    
    # SUB FUNCTION TO PLOT THE DATA
    plottarget <- function(siardata,isox,isoy,a,grps){
    # Plot the target data
    pchseq <- c(1:2,4:20)
        for (j in 1:nrow(siardata$targets)) {
            
            if(siardata$numgroups != 1 & !is.null(grps)){
                if(any(siardata$targets[j,1]==grps)){
                    dx <- siardata$targets[j, isox+a] # - siardata$corrections[isox, 2]
                    dex <- 2 * siardata$corrections[isox, 3]
                    
                    dy <- siardata$targets[j, isoy+a] # - siardata$corrections[isoy, 2]
                    dey <- 2 * siardata$corrections[isoy, 3]
                    
                    siaraddcross(x=dx,y=dy, clr="grey50",upch=pchseq[siardata$targets[j,1]])
                }
            }
            else{
                    dx <- siardata$targets[j, isox+a] # - siardata$corrections[isox, 2]
                    dex <- 2 * siardata$corrections[isox, 3]
                    
                    dy <- siardata$targets[j, isoy+a] #- siardata$corrections[isoy, 2]
                    dey <- 2 * siardata$corrections[isoy, 3]
                    
                    siaraddcross(x=dx,y=dy, clr="grey50")
            }
        }   
    }
    
    # SUB FUNTION THAT HANDLES ALL THE PLOTTING AND CALLS THE OTHER SUB-FUNCTIONS
    fun1 <- function(siardata, siarversion = 0, grp=NULL,panel=NULL,isos=c(1,2)){
        if(!is.null(panel) & is.null(grp)){
            warning(cat("WARNING. grp set to ALL and panel set to a value.\n Overriding your panel selection and setting to panel=NULL.\n In order to plot all groups on seperate panels please call\n grp=1:siardata$numgroups and panel=1 or panel=c(r,c)\n to specify number of rows and columns"))
            panel<-NULL
        }
        # need to select the isotopes here
        if(all(isos==0)){
            isox <- 1
            isoy <- 2
        }
        else{
            isox <- isos[1]
            isoy <- isos[2]
        }
        
        
        # this line takes care of the case when there is only a single group
        a <- 1
        if(siardata$numgroups == 1){a<-0}
        
        # sort out the panels if seperating groups
        if(!is.null(panel)){
            if(prod(panel)<length(grp)){
                panel <- c(ceiling(sqrt(length(grp))))
                panel <- c(max(panel,1),max(ceiling(length(grp)/panel),1))
            }
            
            split.screen(panel)
        }
        else{
            newgraphwindow()
        }
        
        # define the axes limits
        er <- (siardata$sources[, (2*isox)+1]^2 + siardata$corrections[, (2*isox)+1]^2 )^0.5
        xmins <- min(c(siardata$sources[, 2*isox]+ siardata$corrections[, 2*isox] - 3 * er, siardata$targets[, isox+a]))
        xmaxs <- max(c(siardata$sources[, 2*isox]+ siardata$corrections[, 2*isox] + 3 * er, siardata$targets[, isox+a]))
        
        er <- ( siardata$sources[, (2*isoy)+1]^2 + siardata$corrections[, (2*isoy)+1]^2 )^0.5
        ymins <- min(c(siardata$sources[, 2*isoy]+ siardata$corrections[, 2*isoy] - 3 * er, siardata$targets[, isoy+a]))
        ymaxs <- max(c(siardata$sources[, 2*isoy]+ siardata$corrections[, 2*isoy] + 3 * er, siardata$targets[, isoy+a]))

            
        if(is.null(panel)){
            #if (siardata$corrections[1, 1] == 0) 
                plot(1, 1, type = "n", xlim = c(xmins, xmaxs), 
                  ylim = c(ymins, ymaxs), main = siardata$TITLE, 
                  xlab = colnames(siardata$targets)[isox+a], ylab = colnames(siardata$targets)[isoy+a])
            #if (siardata$corrections[1, 1] != 0) 
            #    plot(1, 1, type = "n", xlim = c(xmins, xmaxs), 
            #      ylim = c(ymins, ymaxs), main = siardata$TITLE, 
             #     xlab = siardata$corrections[isox, 1], ylab = siardata$corrections[isoy, 
            #        1])
        }
                    
        
    
                
        for (k in 1:length(grp)){
            
            if(!is.null(panel)){
                screen(k)
            
    
                #if (siardata$corrections[1, 1] == 0) 
                    plot(1, 1, type = "n", xlim = c(xmins, xmaxs), 
                      ylim = c(ymins, ymaxs), main = paste("Group",grp[k]), 
                      xlab = colnames(siardata$targets)[isox+a], ylab = colnames(siardata$targets)[isoy+a])
                #if (siardata$corrections[1, 1] != 0) 
                 #   plot(1, 1, type = "n", xlim = c(xmins, xmaxs), 
                 #     ylim = c(ymins, ymaxs), main = paste("Group",grp[k]), 
                 #     xlab = siardata$corrections[isox, 1], ylab = siardata$corrections[isoy, 
                 #       1])
            }
            
            
            ## HERE          
            # Plot the target data
            if(!is.null(grp)){
                plottarget(siardata,isox,isoy,a,grps=grp[k])
            }
            else{
                plottarget(siardata,isox,isoy,a,grps=grp)
            }
    
            
            
            # Plot the sources
            
            for (i in 1:nrow(siardata$sources)) {
                
                dx <- siardata$sources[i, 2*isox] + siardata$corrections[i, 2*isox]
                dex <- 2 * (siardata$sources[i, (2*isox)+1]^2 + siardata$corrections[i, (2*isox)+1]^2 )^0.5
                
                dy <- siardata$sources[i, 2*isoy] + siardata$corrections[i, 2*isoy]
                dey <- 2 * ( siardata$sources[i, (2*isoy)+1]^2 + siardata$corrections[i, (2*isoy)+1]^2 )^0.5
                
                siaraddcross(x=dx,ex=dex,y=dy,ey=dey,upch=15, clr=i)
            }
        
    
        }
        
        if(!is.null(panel)){
            close.screen(all = TRUE) 
        }
        
        if (siarversion > 0){
                mtext(paste("siar v", siarversion), side = 1, 
                  line = 4, adj = 1, cex = 0.6)
        }
        
        if(siardata$numgroups==0){
            grp<-1
        }
        
        if(is.null(grp)){
            grp<-1
        }
        
        pchseq <- c(1:2,4:20)
        
        # NEED TO CHECK THAT THE LEGEND FUNCTION WORKS WHEN GRP=NULL AND WHEN NUMGROUPS=0
        if (leg == 1){
            datalabs <- NULL
            if(siardata$numgroups==1) {
              datalabs <- "data"
            } else {
              for(k in 1:length(grp)){
                  datalabs <- c(datalabs,as.character(paste("Group",grp[k])))
              }
            }
            legend(locator(1), legend = c(as.character(siardata$sources[, 
                1]), datalabs), lty = c(rep(1, nrow(siardata$sources)),rep(-1,length(grp))),
                pch = c(rep(15, nrow(siardata$sources)), 
                pchseq[grp]), col = c(seq(1, nrow(siardata$sources)), rep("grey50",length(grp))), 
                bty = "n")
        }
        
        if (leg == 2){
            datalabs <- NULL
            if(siardata$numgroups==1) {
              datalabs <- "data"
            } else {
              for(k in 1:length(grp)){
                  datalabs <- c(datalabs,as.character(paste("Group",grp[k])))
              }
            }
            newgraphwindow()
            plot(0,0,"n",xaxt="n",yaxt="n",bty="n")
            legend(0,0, legend = c(as.character(siardata$sources[, 
                1]), datalabs), lty = c(rep(1, nrow(siardata$sources)),rep(-1,length(grp))),
                pch = c(rep(15, nrow(siardata$sources)), 
                pchseq[grp]), col = c(seq(1, nrow(siardata$sources)), rep("grey50",length(grp))), 
                bty = "n")
        }

 

}

    # THE ACTUAL COMMANDS FOR THE OVERSEER FUNCTION ARE HERE
    # CALLS FUN1() MULTIPLE TIMES FOR > 2 ISOTOPES
    if(siardata$numiso>2 & all(isos==0)){
        for (i in 1:(siardata$numiso-1)){
            for (j in (i+1):siardata$numiso){
                fun1(siardata,isos=c(i,j))
            }
        }
    }
    else{
        fun1(siardata, siarversion, grp,panel,isos)
    }
    
    
}
