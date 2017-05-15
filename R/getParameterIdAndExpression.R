getParameterIdAndExpression <- function(parameterXMLNode) {
  if (!any("Id"==XML::xmlSApply(parameterXMLNode, XML::xmlName))) {
    stop("ERROR: Parameter without Id tag")
  }
  if (!any("Expression"==XML::xmlSApply(parameterXMLNode, XML::xmlName))) {
    stop("ERROR: Parameter without Expression tag")
  }
  #   # might contain a Description tag?
  #   if (xmlSize(parameterXMLNode)!=2) {
  #     stop("ERROR: Parameter must contain exactly two tags (Id and Expression)")
  #   }
  param <- c(Id=XML::xmlValue(parameterXMLNode[["Id"]]),Expression=XML::xmlValue(parameterXMLNode[["Expression"]]))
  if (sub("[[:space:]]+$", "", param["Id"]!=param["Id"])) {
    warning(paste("White space detected in Parameter Id ",param["Id"]," removing spaces"))
    param["Id"]=sub("[[:space:]]+$", "", param["Id"])
  }
  return(param)
}

