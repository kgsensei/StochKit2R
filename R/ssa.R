#'@title Gillespie Stochastic Simulation Algorithm
#'
#'@description
#'\code{ssa} Runs an ensemble of SSA \code{realizations} and writes output
#'data to \code{outputDir}
#'
#'@param modelFile Character string with path to StochKit2 .xml model file
#'@param time Simulation time of each realization
#'@param realizations Number of realizations
#'@param intervals Number of output intervals. Default 0 outputs at end time only. 1=keep data at start and end time, 2=keep data at start, middle, and end times, etc. Note data is stored at (intervals+1) equally spaced time points.
#'@param noStats Do not keep means and variances data. Default \code{FALSE} creates stats directory in output directory
#'@param keepTrajectories Keep trajectory data. Creates trajectories directory in output directory
#'@param keepHistograms Keep histogram data. Creates histograms directory in output directory
#'@param bins Number of histogram bins
#'@param outputDir Character string with path to output directory. By default (=NULL) and will not write file. If output directory does not exist, it will be created. If output directory already exists, use \code{force=TRUE} to overwrite
#'@param force Force overwriting of existing data
#'@param seed Seed the random number generator. By default the seed is determined by the R random number generator, so the seed can also be set by calling \code{set.seed} in R immediately before calling \code{ssa}
#'@param p Override default and specify the number of processes (threads) to use. By default (=0), the number of processes will be determined automatically (recommended). Ignored on systems without OpenMP support.
#'@return List of data frames with means, variances, and a list of trajectory vectors
#'@examples
#'\dontrun{
#'#'#example using included dimer_decay.xml file
#'#output written to directory ex_out (created in current working directory)
#'#run 100 simulations for 10 time units, keeping output at 20 time intervals
#'#store model file name in a variable first
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'out <- ssa(model,"ex_out",10,100,20,force=TRUE)
#'
#'#more typical example where model file is stored elsewhere
#'#(must be valid path to existing .xml StochKit2 model file)
#'#store output in dimer_decay_output, overwrite existing data
#'#and keep trajectory data.
#'out <- ssa("Desktop/dimer_decay.xml",
#'    "Desktop/dimer_decay_output",10,100,20,keepTrajectories=T,force=T)
#'}
ssa <- function(modelFile,time,realizations,intervals=0,noStats=FALSE,keepTrajectories=FALSE,keepHistograms=FALSE,bins=32,outputDir=NULL,force=FALSE,seed=NULL,p=0) {
  # can set seed in R with set.seed()
  
  #checks on modelFile  
  if (!file.exists(modelFile)) {
    stop("ERROR: Model file does not exist")
  }
  #simple checks
  if (time<=0) {
    stop("ERROR: simulation time must be positive")
  }
  if (realizations<=0) {
    stop("ERROR: realizations must be positive")
  }
  if (intervals<0) {
    stop("ERROR: intervals must not be negative")
  }
  if (bins<=0) {
    stop("ERROR: number of bins must be positive")
  }
  if (p<0) {
    stop("ERROR: number of processes must not be negative (when p=0, processes will be determined automatically)")
  }
  
  
  #must keep some output
  if (noStats & !keepTrajectories & !keepHistograms) {
    stop("ERROR: must output at least one of stats, trajectories, or histograms.")
  }
  
  
  if(!is.null(outputDir)){
    #expand path
    outputDir <- tryCatch(suppressWarnings(normalizePath(outputDir)), error = function(e) {stop("Invalid or missing outputDir output directory path, terminating StochKit2R")}, finally = NULL)
    #remove tailing slashes or backslashes
    #because file.exists returns false if directory ends in slash or backslash
    outputDir <- gsub("//*$","",outputDir)
    outputDir <- gsub("\\\\*$","",outputDir)

    createOutputDirs(outputDir,noStats,keepTrajectories,keepHistograms,force)
  }
  if(is.null(outputDir)){
    outputDir=""
  }
  
  if (is.null(seed)) {
    seed=floor(runif(1,-.Machine$integer.max,.Machine$integer.max))
  }
  
  # process the xml model file, put in R Lists
  StochKit2Rmodel <- buildStochKit2Rmodel(modelFile)
  
  # everything is set up, ready to run the simulation
  ssaStochKit2RInterface(StochKit2Rmodel,time,realizations,intervals,!noStats,keepTrajectories,keepHistograms,bins,outputDir,seed,p)
}