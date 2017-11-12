## ----, eval=FALSE--------------------------------------------------------
#  citation(package="StochKit2R")

## ----, eval=FALSE--------------------------------------------------------
#  install.packages("StochKit2R")

## ----, eval=FALSE--------------------------------------------------------
#  library(StochKit2R)

## ----, eval=FALSE--------------------------------------------------------
#  install.packages("<path to>/StochKit2R_<version>.tar.gz",repos=NULL,type="source")

## ----eval=FALSE----------------------------------------------------------
#  # use the dimer_decay.xml model that is included with StochKit2R
#  # to make things easier to read, store the path in a variable
#  model <- system.file("dimer_decay.xml",package="StochKit2R")
#  # now, run a simulation, you can save the output locally and save it to a selected output:
#  out <- ssa(model,"ex_out",10,100,20)

## ----, eval=FALSE--------------------------------------------------------
#  out <- ssa("~/Desktop/dimer_decay.xml","~/Desktop/dimer_decay_output",10,100,20)

## ----, eval=FALSE--------------------------------------------------------
#  out <- ssa("~/Desktop/dimer_decay.xml","~/Desktop/dimer_decay_output",10,100,20,force=T)

## ----, eval=FALSE--------------------------------------------------------
#  out <- ssaSingle(system.file("dimer_decay.xml",package="StochKit2R"),"single_output.txt",0,10)

## ----, fig.width=6, fig.height=4, fig.align="center"---------------------
library(StochKit2R)
#example using included dimer_decay.xml file
model <- system.file("dimer_decay.xml",package="StochKit2R")
#output list saved to out
out <- ssa(model,"ex_out",10,100,20,force=TRUE)
#plot the data for species 2 and 3 (all of them in the dimer decay model)
plotStats(out,c(2,3))

## ----, eval=FALSE--------------------------------------------------------
#  model <- system.file("dimer_decay.xml",package="StochKit2R")
#  file.copy(model,"~/Desktop/dimer_decay.xml")

## ----, comment=""--------------------------------------------------------
model <- system.file("dimer_decay.xml",package="StochKit2R")
XML::xmlParse(model)

