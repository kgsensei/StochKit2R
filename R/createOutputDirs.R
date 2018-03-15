#'@title Create Output Directories for ssa and tauLeaping
#'
#'@description
#'\code{createOutputDirs} Called by writeOutput create output directory
#'
#'@param outputDir Character string with path to output directory.
#'@param subDirs Character vector
#'@param force delete and overwrite existing output directories, if they exists
createOutputDirs <- function(outputDir,subdirs,force=FALSE) {
  if (file.exists(outputDir)) {
    if (!force) {
      stop("ERROR: Output directory already exists. Delete existing directory, specify a unique directory name, or run with force=TRUE to overwrite")
    }
    else {
      # delete subdirectories
      #lapply(dir(outputDir),function(f) unlink(file.path(outputDir,f),recursive=T))
      fail = unlink(file.path(outputDir,dir(outputDir)),recursive=T)==1
      if (fail) {
        stop("ERROR: unable to delete existing output directory (or subdirectory)")
      }
    }
  } else {
  # create outputDir
    fail = !dir.create(file.path(outputDir),recursive=TRUE)
    if (fail) {
      stop("ERROR: unable to create output directory")
    }
  }
  
  # create subdirs
  fail = !all(sapply(subdirs,function(d) dir.create(file.path(outputDir,d),recursive=T)))
  if (fail) {
    stop("ERROR: unable to create all subdirectories")
  }
  
}
  
  