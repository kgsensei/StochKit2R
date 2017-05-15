#note: InitialPopulation should be treated as a string since it
#might contain parameters that need to be evaluated
getSpeciesIdAndInitialPopulation <- function(speciesXMLNode) {
  if (!any("Id"==XML::xmlSApply(speciesXMLNode, XML::xmlName))) {
    stop("ERROR: Species without Id tag")
  }
  if (!any("InitialPopulation"==XML::xmlSApply(speciesXMLNode, XML::xmlName))) {
    stop("ERROR: Parameter without InitialPopulation tag")
  }
  spec <- c(Id=XML::xmlValue(speciesXMLNode[["Id"]]),InitialPopulation=XML::xmlValue(speciesXMLNode[["InitialPopulation"]]))
  if (sub("[[:space:]]+$", "", spec["Id"])!=spec["Id"]) {
    warning(paste("White space detected in Species Id ",spec["Id"],", removing spaces"))
    spec["Id"]=sub("[[:space:]]+$", "", spec["Id"])
  }
  return(spec)
}

