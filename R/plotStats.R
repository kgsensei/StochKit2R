#'@title Plot StochKit2R simulation statistics output data
#'
#'@description
#'\code{plotStats} Plots means and mean +/- one standard deviation of populations specified in \code{indices} in the StochKit2R stats output directory \code{statsDirectory}
#'
#'@param statsDirectory Character string with path to StochKit2 stats output directory
#'@param indices The species indices that will be plotted. The first species is index 1
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'ssa(model,"ex_out",10,100,20,force=TRUE)
#'#plot the data for species 1,2 and 3 (all of them in the dimer decay model)
#'plotStats("ex_out/stats",c(1,2,3))
#'}
plotStats <- function(statsDirectory,indices) {
  
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
 
  # get file names
  fnameMean = paste(statsDirectory,'/means.txt',sep='');
  fnameVar =  paste(statsDirectory,'/variances.txt',sep='');
  
  #read the first line of the means file
  #and check for headers (labels)
  line1 <- strsplit(readLines(fnameMean,n=1),split="\t")[[1]]
  if (line1[1]=="time") {
    hasLabels=TRUE
  } else {
    hasLabels=FALSE
  }
  
  #get the means data
  meansData <- read.table(fnameMean,header=hasLabels)[,c(1,indices)]
  #give (slightly) more meaningful labels if none
  if (!hasLabels) {
    names(meansData) <- c("time",names(meansData)[1:(length(meansData)-1)])
  }
  #get the variances data
  variancesData <- read.table(fnameVar,header=hasLabels)[,c(1,indices)]
  #give (slightly) more meaningful labels if none
  if (!hasLabels) {
    names(variancesData) <- names(meansData)
  }
  
  #put data into plottable form
  value=NULL#only to appease R CMD CHECK for CRAN
  variable=NULL#only to appease R CMD CHECK for CRAN
  meansDataMelted <- melt(meansData,id="time")
  p <- qplot(time,value,data=meansDataMelted,colour=variable) + geom_line() + ylab("population") + ggtitle("stats plot")
  # create mean+stdev data
  plusStdevData <- meansData
  plusStdevData[,2:ncol(plusStdevData)] <- plusStdevData[,2:ncol(plusStdevData)] + sqrt(variancesData[,2:ncol(variancesData)])
  plusStdevDataMelted <- melt(plusStdevData,id="time")
  p2 <- p + geom_line(data=plusStdevDataMelted,linetype="dashed")
  # create mean-stdev data
  minusStdevData <- meansData
  minusStdevData[,2:ncol(minusStdevData)] <- minusStdevData[,2:ncol(minusStdevData)] - sqrt(variancesData[,2:ncol(variancesData)])
  minusStdevDataMelted <- melt(minusStdevData,id="time")
  p3 <- p2 + geom_line(data=minusStdevDataMelted,linetype="dashed")
  return(p3)
}