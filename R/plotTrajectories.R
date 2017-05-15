#'@title Plot StochKit2R simulation trajectory output data
#'
#'@description
#'\code{plotTrajectories} Plots trajectories numbered \code{startn} to \code{endn} of populations specified in \code{indices} in the StochKit2R trajectories output directory \code{trajectoriesDirectory}
#'
#'@param trajectoriesDirectory Character string with path to StochKit2 trajectories output directory
#'@param startn First trajectory to print (note: trajectories start at index 0)
#'@param endn Last trajectory to print
#'@param indices Vector of the species indices that will be plotted. The first species is index 1
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'ssa(model,"ex_out",10,100,20,keepTrajectories=TRUE,force=TRUE)
#'#plot the data for species 1,2 and 3 for trajectories 2,3,4 and 5
#'plotTrajectories("ex_out/trajectories",2,5,c(1,2,3))
#'}
plotTrajectories <- function(trajectoriesDirectory,startn,endn,indices) {
  
  if ((round(startn)!=startn) | (round(endn) != endn)) {
    stop('Need integer numbers for input file indexing')
  }
  if (startn > endn) {
    stop('start file index is larger than the end file index')
  }
  if (length(indices)==0) {
    stop('indices vector must not be empty')
  }
  
  # Check to make sure that the indices are integers
  if (sum((round(indices)==indices)) != length(indices)) {
    stop('Indexes must be integers')
  }
  
  # Check to make sure that the indices are greater than 0
  if (min(indices)<=0) {
    stop('Indexes must be positive integers')
  }
  
  # Check to eliminate duplicate indices
  indices2 = unique(indices);
  if (length(indices2)!=length(indices)) {
    warning(paste("removing ",length(indices) - length(indices2), " duplicate indexes."))
    indices=indices2;
  }
  
  indices = indices+1; # add 1 to account for the time column
  indices = sort(indices);
  len_indices =length(indices);
 
  #read the first line of the first file
  #and check for headers (labels)  
  fileName <- paste(trajectoriesDirectory,"/trajectory",startn,".txt",sep="")
  line1 <- strsplit(readLines(fileName,n=1),split="\t")[[1]]
  if (line1[1]=="time") {
    hasLabels=TRUE
  }
  else {
    hasLabels=FALSE
  }
  
  #plot data...
  fileName <- paste(trajectoriesDirectory,"/trajectory",startn,".txt",sep="")
  trajData <- read.table(fileName,header=hasLabels)[,c(1,indices)]
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
  if (startn < endn) {
    for (i in (startn+1):endn) {
      fileName <- paste(trajectoriesDirectory,"/trajectory",i,".txt",sep="")
      trajData <- read.table(fileName,header=hasLabels)[,c(1,indices)]
      #give (slightly) more meaningful labels if none
      if (!hasLabels) {
        names(trajData) <- c("time",names(trajData)[1:(length(trajData)-1)])
      }
      
      trajDataMelted2 <- melt(trajData,id="time")
      p <- p + geom_line(data=trajDataMelted2)
    }
  }
  return(p)
}
