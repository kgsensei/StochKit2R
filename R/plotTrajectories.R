#'@title Plot StochKit2R simulation trajectory output data
#'
#'@description
#'\code{plotTrajectories} Plots trajectories \code{outputIndexes} of populations specified in \code{speciesIndexes} in the StochKit2R trajectories in \code{trajectoriesData}
#'
#'@param trajectoriesData trajs list from ssa (or tauLeaping) returned output or character string containing the path to the trajectories output directory
#'@param outputIndexes Integer vector containing indexes for which trials will be plotted
#'@param speciesIndexes Integer vector of the species indices that will be plotted. The first species is index 1
#'@param file set to TRUE if trajectoriesData is a path to the trajectories output directory (rather than returned output data)
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,time=10,realizations=100,intervals=20,
#'           keepTrajectories=TRUE,outputDir="ex_out",force=TRUE)
#'#plot the data for species 1,2 and 3 for trajectories 2,3,4 and 5 from file
#'plotTrajectories("ex_out/trajectories",outputIndexes=2:3,speciesIndexes=1:3,file=TRUE)
#'#same plot from output
#'plotTrajectories(out$trajs,outputIndexes=2:3,speciesIndexes=1:3)
#'}
plotTrajectories <- function(trajectoriesData,outputIndexes,speciesIndexes,file=F) {

  #check to make sure indices are provided
  if (length(outputIndexes)==0 | length(speciesIndexes)==0) {
    stop('index vectors must not be empty')
  }
  
  # Check to make sure that the indices are integers
  if (sum(round(outputIndexes)!=outputIndexes)!=0 | sum(round(speciesIndexes)!=speciesIndexes)!=0) {
    stop('indices must be integers')
  }
  
  # Check to make sure that the indices are greater than 0
  if (sum(outputIndexes<=0)!=0 | sum(speciesIndexes<=0)!=0) {
     stop('all indices must be greater than 0')
  }
  
  # Check to eliminate duplicate indices
  indices2 = unique(outputIndexes);
  if (length(indices2)!=length(outputIndexes)) {
    warning(paste("removing ",length(outputIndexes) - length(indices2), " duplicate indexes from outputIndexes."))
    outputIndexes=indices2;
  }
  
  indices2 = unique(speciesIndexes);
  if (length(indices2)!=length(speciesIndexes)) {
    warning(paste("removing ",length(speciesIndexes) - length(indices2), " duplicate indexes from speciesIndexes."))
    speciesIndexes=indices2;
  }

  speciesIndexes = speciesIndexes+1; # add 1 to account for the time column
  speciesIndexes = sort(speciesIndexes);
  
  if(!file){
    #check for headers (labels)

    colLabels <- names(trajectoriesData[[1]])
    if (colLabels[1]=="time") {
      hasLabels=TRUE
    }
    else {
      hasLabels=FALSE
    }

  #plot data...
    
    trajData <- trajectoriesData[[1]][,c(1,speciesIndexes)]
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
    if (length(outputIndexes>1)) {
      for (i in outputIndexes[2:length(outputIndexes)]) {
        trajData <- trajectoriesData[[i]][,c(1,speciesIndexes)]
        #give (slightly) more meaningful labels if none
        if (!hasLabels) {
          names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
        }
        
        trajDataMelted2 <- melt(trajData,id="time")
        p <- p + geom_line(data=trajDataMelted2)
      }
     }  

  }

  else{
    #read the first line of the first file
    #and check for headers (labels)  
    fileName <- paste(trajectoriesData,"/trajectory",outputIndexes[1]-1,".txt",sep="")
    line1 <- strsplit(readLines(fileName,n=1),split="\t")[[1]]
    if (line1[1]=="time") {
      hasLabels=TRUE
    }
    else {
      hasLabels=FALSE
    }

  #plot data...
    fileName <- paste(trajectoriesData,"/trajectory",outputIndexes[1]-1,".txt",sep="")
    trajData <- read.table(fileName,header=hasLabels)[,c(1,speciesIndexes)]
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
    if (length(outputIndexes)>1) {
      for (i in outputIndexes[2:length(outputIndexes)-1]) {
        fileName <- paste(trajectoriesData,"/trajectory",outputIndexes[i]-1,".txt",sep="")
        trajData <- read.table(fileName,header=hasLabels)[,c(1,speciesIndexes)]
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
