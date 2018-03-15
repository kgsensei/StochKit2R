## ---- eval=FALSE---------------------------------------------------------
#  citation(package="StochKit2R")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("StochKit2R")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages(file.choose(),repos=NULL)

## ---- eval=FALSE---------------------------------------------------------
#  library(StochKit2R)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools") #run once on your system
#  devtools::install_github("StochKit2R/StochKit2R")

## ---- eval=FALSE---------------------------------------------------------
#  library(StochKit2R)

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("devtools") #run once on your system
#  devtools::install_github("StochKit2R/StochKit2R",ref="serial")

## ---- eval=FALSE---------------------------------------------------------
#  install.packages("<path to>/StochKit2R_<version>.tar.gz",repos=NULL,type="source")

## ----eval=FALSE----------------------------------------------------------
#  # use the dimer_decay.xml model that is included with StochKit2R
#  # to make things easier to read, store the path in a variable
#  model = system.file("dimer_decay.xml",package="StochKit2R")
#  # now, run a simulation:
#  out = ssa(modelFile=model,time=10,realizations=100,intervals=20)

## ---- eval=FALSE---------------------------------------------------------
#  out = ssa("~/Desktop/dimer_decay.xml",10,100,20)

## ---- eval=FALSE---------------------------------------------------------
#  out = ssaSingle(system.file("dimer_decay.xml",package="StochKit2R"),10)

## ---- eval=FALSE---------------------------------------------------------
#  out = ode(system.file("dimer_decay.xml",package="StochKit2R"),time=10,intervals=20)

## ---- fig.width=6, fig.height=4, fig.align="center"----------------------
library(StochKit2R)
#example using included dimer_decay.xml file
model = system.file("dimer_decay.xml",package="StochKit2R")
##output value stored locally for use and written to ex_out file in working directory
out = ssa(modelFile=model,time=10,realizations=100,interval=20)
#plot the data for species 2 and 3 (indexes; could also use names)
plotStats(out,c(2,3))

## ---- eval=FALSE---------------------------------------------------------
#  model = system.file("dimer_decay.xml",package="StochKit2R")
#  file.copy(model,"~/Desktop/dimer_decay.xml")

## ---- comment=""---------------------------------------------------------
model = system.file("dimer_decay.xml",package="StochKit2R")
writeLines(readLines(model))

