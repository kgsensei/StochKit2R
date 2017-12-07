#'@title Gillespie Stochastic Simulation Algorithm single trajectory
#'
#'@description
#'\code{ssaSingle} Run a single SSA trajectory. Stores a row of data for every reaction event.
#'
#'@param modelFile Character string with path to StochKit2 .xml model file
#'@param startTime Simulation start time
#'@param endTime Simulation end time
#'@param outputFile Character string with path to output file. By default (=NULL) no data is written to file.
#'@param seed Seed the random number generator. By default the seed is determined by the R random number generator, so the seed can also be set by calling \code{set.seed} in R immediately before calling \code{ssaSingle}
#'@return Data frame containing time and species population for every reaction event in the simulation.
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'#output written to file single_output.txt (created in current working directory)
#'out <- ssaSingle(system.file("dimer_decay.xml",package="StochKit2R"),startTime=0,endTime=10)
#'}
ssaSingle <- function(modelFile,startTime,endTime,outputFile=NULL,seed=NULL) {
  # can set seed in R with set.seed()
  
  #checks on modelFile  
  if (!file.exists(modelFile)) {
    stop("ERROR: Model file does not exist")
  }
  #simple checks
  if (startTime>=endTime) {
    stop("startTime must be before endTime")
  }
  
  #checks on outputDirand expand path
  if(!is.null(outputFile)){
    outputFile <- tryCatch(suppressWarnings(normalizePath(outputFile)), error = function(e) {stop("Invalid or missing outputFile path, terminating StochKit2R")}, finally = NULL)
    if (file.exists(outputFile)) {
        stop("Output file already exists. Delete existing file or specify a unique file name")
    }
  }
    
  if (is.null(seed)) {
    seed=floor(runif(1,-.Machine$integer.max,.Machine$integer.max))
  }
  
  #changes null ourputDir to blank string because c++ strings are not nullable
  if(is.null(outputFile)){
    outputFile=""
  }

  # process the xml model file, put in R Lists
  StochKit2Rmodel <- buildStochKit2Rmodel(modelFile)
  
  # everything is set up, ready to run the simulation
  ssaSingleStochKit2RInterface(StochKit2Rmodel,startTime,endTime,outputFile,seed)
}