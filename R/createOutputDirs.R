#'@title Create Output Directories for ssa and tauLeaping
#'
#'@description
#'\code{createOutputDirs} Called to ssa and tauLeaping to create output directory
#'
#'@param outputDir Character string with path to output directory.
#'@param noStats true if not keeping statistics data
#'@param keepTrajectories true if keeping trajectory data
#'@param keepHistograms true if keeping histogram data
#'@param overwrite and delete existing output directories, if they exists
createOutputDirs <- function(outputDir,noStats,keepTrajectories,keepHistograms,force) {
  if (file.exists(outputDir)) {
    if (!force) {
      stop("ERROR: Output directory already exists. Delete existing directory, specify a unique directory name, or run with force=TRUE to overwrite")
    }
    else {
      #delete stats, trajectories, histograms directories, if they exist
      unlink(file.path(outputDir,"stats"),recursive=TRUE)
      unlink(file.path(outputDir,"trajectories"),recursive=TRUE)
      unlink(file.path(outputDir,"histograms"),recursive=TRUE)
    }
  }
  else {
    dir.create(file.path(outputDir),recursive=TRUE)
    #mkdirs(outputDir)
  }
  
  
  #create stats, trajectories, histograms output directories, if needed
  if (!noStats) {
    dir.create(file.path(outputDir,"stats"))
  }
  if (keepTrajectories) {
    dir.create(file.path(outputDir,"trajectories"))
  }
  if (keepHistograms) {
    dir.create(file.path(outputDir,"histograms"))
  }
  
}
  
  