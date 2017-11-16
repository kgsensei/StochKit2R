#'@title Plot StochKit2R simulation trajectory output data
#'
#'@description
#'\code{plotTrajectories} Plots trajectories \code{outputIndex} of populations specified in \code{speciesIndex} in the StochKit2R trajectories in \code{trajectoriesData}
#'
#'@param trajectoriesData String containing file (must end in /trajectories) name or ssa output, depending on the value of the file parameter
#'@param outputIndex Integer vector containing indexes for which trials will be plotted
#'@param speciesIndex Integer Vector of the species indices that will be plotted. The first species is index 1
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,"ex_out",10,100,20,keepTrajectories=TRUE,force=TRUE)
#'#plot the data for species 1,2 and 3 for trajectories 2,3,4 and 5
#'plotTrajectories("ex_out/trajectories",2:3,1:3,file=T)
#'}
plotTrajectories <- function(trajectoriesData,outputIndex,speciesIndex,file=F) {

  
  if (length(outputIndex)==0 | length(speciesIndex)==0) {
    stop('index vectors must not be empty')
  }
  
  # Check to make sure that the indices are integers
  if (sum(round(outputIndex)!=outputIndex)!=0 | sum(round(speciesIndex)!=speciesIndex)!=0) {
    stop('indices must be integers')
  }
  
  # Check to make sure that the indices are greater than 0
  if (sum(outputIndex<=0)!=0 | sum(speciesIndex<=0)!=0) {
     stop('all indices must be greater than 0')
  }
  
  # Check to eliminate duplicate indices
  indices2 = unique(outputIndex);
  if (length(indices2)!=length(outputIndex)) {
    warning(paste("removing ",length(outputIndex) - length(indices2), " duplicate indexes from outputIndex."))
    outputIndex=indices2;
  }
  
  indices2 = unique(speciesIndex);
  if (length(indices2)!=length(speciesIndex)) {
    warning(paste("removing ",length(speciesIndex) - length(indices2), " duplicate indexes from speciesIndex."))
    speciesIndex=indices2;
  }

  speciesIndex = speciesIndex+1; # add 1 to account for the time column
  speciesIndex = sort(speciesIndex);
  
  if(!file){
    #check for headers (labels)
    trajs=trajectoriesData$trajs

    colLabels <- names(trajs[[1]])
    if (colLabels[1]=="time") {
      hasLabels=TRUE
    }
    else {
      hasLabels=FALSE
    }

  #plot data...
    
    trajData <- trajs[[1]][,c(1,speciesIndex)]
      #give (slightly) more meaningful labels if none
    if (!hasLabels) {
      names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
    }

    #put data into plottable form
    value=NULL#only to appease R CMD CHECK for CRAN
    variable=NULL#only to appease R CMD CHECK for CRAN
    trajDataMelted <- melt(trajData,id="time")
    p <- qplot(time,value,data=trajDataMelted,colour=variable,geom="line") + ylab("population") + ggtitle("trajectory plot")
    
    # plot the rest of the trajectories!
    if (length(outputIndex>1)) {
      for (i in outputIndex[2:length(outputIndex)]) {
        trajData <- trajs[[i]][,c(1,speciesIndex)]
        #give (slightly) more meaningful labels if none
        if (!hasLabels) {
          names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
        }
        
        trajDataMelted2 <- melt(trajData,id="time")
        p <- p + geom_line(data=trajDataMelted2)
      }
     }  

  }



  if(file){
    #read the first line of the first file
    #and check for headers (labels)  
    fileName <- paste(trajectoriesData,"/trajectory",outputIndex[1]-1,".txt",sep="")
    line1 <- strsplit(readLines(fileName,n=1),split="\t")[[1]]
    if (line1[1]=="time") {
      hasLabels=TRUE
    }
    else {
      hasLabels=FALSE
    }

  #plot data...
    fileName <- paste(trajectoriesData,"/trajectory",outputIndex[1]-1,".txt",sep="")
    trajData <- read.table(fileName,header=hasLabels)[,c(1,speciesIndex)]
      #give (slightly) more meaningful labels if none
      if (!hasLabels) {
        names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
      }

    #put data into plottable form
    value=NULL#only to appease R CMD CHECK for CRAN
    variable=NULL#only to appease R CMD CHECK for CRAN
    trajDataMelted <- melt(trajData,id="time")
    p <- qplot(time,value,data=trajDataMelted,colour=variable,geom="line") + ylab("population") + ggtitle("trajectory plot")
    
    # plot the rest of the trajectories!
    if (length(outputIndex)>1) {
      for (i in outputIndex[2:length(outputIndex)-1]) {
        fileName <- paste(trajectoriesData,"/trajectory",outputIndex[i]-1,".txt",sep="")
        trajData <- read.table(fileName,header=hasLabels)[,c(1,speciesIndex)]
        #give (slightly) more meaningful labels if none
        if (!hasLabels) {
          names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
        }
        
        trajDataMelted2 <- melt(trajData,id="time")
        p <- p + geom_line(data=trajDataMelted2)
      }
    }
  }
return(p)
}
