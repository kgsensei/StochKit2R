#'@title Plot StochKit2R simulation statistics output data
#'
#'@description
#'\code{plotStats} Plots means and mean +/- one standard deviation of populations specified in \code{indices}
#'
#'@param statsData stats element from ssa (or tauLeaping) output (a list containing the means and vars (variances) data frames) or a character string to the path to the output stats directory.
#'@param species A character vector of species names or a numeric vector of species indexes of the species that will be plotted. For numeric indexes, the first species is index 1. By default =NULL, all species are plotted.
#'@param file set to TRUE if statsData is a path to the stats output directory (rather than returned output data)
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'out <- ssa(model,time=10,realizations=100,intervals=20,outputDir="ex_out",force=TRUE)
#'#plot the data for all species
#'plotStats(out$stats)
#'#plot species "S1" and "S3" using data from file
#'plotStats("ex_out/stats",c(1,2,3),file=TRUE)
#'#same plot from returned output
#'}
plotStats <- function(statsData,species=NULL,file=FALSE) {

  # if user provided the entire output object, help them
  # check if it is data vs file name
  if (class(statsData)=="list") { 
    if (is.null(statsData$means)) {
      if (is.null(statsData$stats$means)) {
        stop("ERROR: invalid statsData")
      } else {
        statsData = statsData$stats
        message("Using statsData$stats")
      }
    }
  }
  
  if (!is.null(species)) {
    if (! (class(species)=="character" || class(species)=="numeric" || class(species)=="integer")) {
      stop("Species must be a character vector of species names or a numeric vector")
    }
    
    if (class(species)=="numeric" || class(species)=="integer") {
      # Check to make sure that the indices are integers
      if (sum((round(species)==species)) != length(species)) {
        stop('Species indexes must be integers')
      }
      
      # Check to make sure that the indices are greater than 0
      if (min(indices)<=0) {
        stop('Species indexes must be positive integers')
      }
    }
    
    # Check to eliminate duplicate species
    species2 = unique(species);
    if (length(species2)!=length(species)) {
      warning(paste("removing ",length(species) - length(species2), " duplicate species."))
      species=species2;
    }
    
    #if numeric, add 1 (if character, handle it later)
    if (class(species)=="numeric" || class(species)=="integer") {
      indices = species+1; # add 1 to account for the time column
      indices = sort(indices);
      len_indices =length(indices);
    }
  }
  
  # by here, species might be NULL
if(!file){
  if (names(statsData$means[1])=="time") {
    hasLabels=TRUE
  } else {
    hasLabels=FALSE
  }
    
    if (class(species)=="character" && !hasLabels) {
      stop("ERROR: data does not have labels, so species names cannot be used. Use species indexes instead.")
    }
    #convert species names to indices
    if (class(species)=="character") {
      indices = sapply(species,function(s) which(s==names(statsData$means)))
      # if one or more names is not found
      if (class(indices)=="list") {
        bad_names = paste(species[sapply(indices,function(i) length(i)==0)],collapse=",")
        stop(paste("ERROR: invalid species name(s) (",bad_names,") entered.",sep=""))
      }
    }
  
  if (is.null(species)) {
    indices=1:(ncol(statsData$means)-1)+1
  }
    #get the means data
    meansData <- statsData$means[,c(1,indices)]
    #give (slightly) more meaningful labels if none
    if (!(names(statsData$means[1])=="time")) {
      names(meansData) <- c("time",names(meansData)[1:(length(meansData)-1)])
    }
    #get the variances data
    variancesData <- statsData$vars[,c(1,indices)]
    #give (slightly) more meaningful labels if none
    if (!(names(statsData$vars[1])=="time")) {
      names(variancesData) <- names(meansData)
    }
}

else{
    # get file names
    fnameMean = paste(statsData,'/means.txt',sep='')
    fnameVar =  paste(statsData,'/variances.txt',sep='')

    if (!(file.exists(fnameMean) && file.exists(fnameVar))) {
      # means and variances files do not exist
      # see if the user passed the outer directory name
      fnameMean = paste(statsData,'/stats/means.txt',sep='')
      fnameVar =  paste(statsData,'/stats/variances.txt',sep='')
      if (!(file.exists(fnameMean) && file.exists(fnameVar))) {
        stop(paste("ERROR: means and/or variances file not found in directory",statsData))
      } else {
        message("Using statsData/stats")
      }
    }
    #read the first line of the means file
    #and check for headers (labels)
    line1 <- strsplit(readLines(fnameMean,n=1),split="\t")[[1]]
    if (line1[1]=="time") {
      hasLabels=TRUE
      } else {
        hasLabels=FALSE
      }

    if (class(species)=="character" && !hasLabels) {
      stop("ERROR: data does not have labels, so species names cannot be used. Use species indexes instead.")
    }
    #convert species names to indices
    if (class(species)=="character") {
      indices = sapply(species,function(s) which(s==line1))
      # if one or more names is not found
      if (class(indices)=="list") {
        bad_names = paste(species[sapply(indices,function(i) length(i)==0)],collapse=",")
        stop(paste("ERROR: invalid species name(s) (",bad_names,") entered.",sep=""))
      }
    }
    
    if (is.null(species)) {
      indices=1:(length(line1)-1)+1
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