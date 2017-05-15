#'@title Plot StochKit2R simulation histogram output data
#'
#'@description
#'\code{plotHistogram} Plots histogram of data stored in StochKit2R histogram output file \code{histogramFileName}. IMPORTANT: histogram file names have format hist_<species index>_<time point>.dat, where species index STARTS AT 0 (not 1!)
#'
#'@param histogramFileName Character string with path to StochKit2 histogram output file. IMPORTANT: histogram file names have format hist_<species index>_<time point>.dat, where species index STARTS AT 0 (not 1!)
#'@return The ggplot object
#'@examples
#'\dontrun{
#'#example using included dimer_decay.xml file
#'model <- system.file("dimer_decay.xml",package="StochKit2R")
#'#output written to ex_out directory (created in current working directory)
#'ssa(model,"ex_out",10,100,20,keepHistograms=TRUE,force=TRUE)
#'#plot the histogram for species 2 at time point 5
#'#IMPORTANT: histogram file names have format:
#'#hist_<species index>_<time point>.dat, where species index STARTS AT 0
#'plotHistogram("ex_out/histograms/hist_1_5.dat")
#'}
plotHistogram <- function(histogramFileName) {
  
  # read in the lines 
  lines <- strsplit(readLines(histogramFileName),split="\t")
  sID <- lines[[1]][1] # species ID
  time = as.numeric(lines[[1]][2]) # time
  sInd = as.numeric(lines[[1]][2]) # species index
  tInd = as.numeric(lines[[1]][3]) # time index

  lb = as.numeric(lines[[2]][1]) # lower bound
  ub = as.numeric(lines[[2]][2]) # upper bound
  width = as.numeric(lines[[2]][3]) # width
  bsize = as.numeric(lines[[2]][4]) # number of bins

  # read in the third line (bin counts)
  data = as.numeric(lines[[3]])
  # error checking
  if (bsize!=length(data)) {
    stop('retrieved size of bins are not equal to the number of bins in the data line')
  }

  dataRange <- seq(lb,ub,by=width)
  bincenters <- (dataRange[-1] + dataRange[-length(dataRange)])/2
  centers=NULL#only to appease R CMD CHECK for CRAN
  counts=NULL#only to appease R CMD CHECK for CRAN
  df <- data.frame(centers=bincenters,counts=data)
  
  ggplot(data=df,aes(x=centers,y=counts)) + geom_bar(stat="identity",width=width,colour="dark blue",fill="dark blue") + ylab("Bin Counts") + xlab("Bin Centers") + ggtitle(paste("Species:",sID," Time:",time))
}
