# parse data from the SpeciesReference data
# function works on Reactants or Products data
# reactantOrProduct is string "Reactant" or "Product", used for error messages
getReactantOrProduct <- function(reactantsOrProductXMLNode, reactantOrProduct) {
  tmpId=XML::xmlGetAttr(reactantsOrProductXMLNode,"id")
  tmpStoich=XML::xmlGetAttr(reactantsOrProductXMLNode,"stoichiometry")
  # check for valid id
  if (sub("[[:space:]]+$", "", tmpId)!=tmpId) {
    warning(paste("White space detected in ",reactantOrProduct," id ",tmpId,", removing spaces",sep=""))
    tmpId=sub("[[:space:]]+$", "", tmpId)
  }
  # check for valid stoichiometry (int)
  if (is.na(as.numeric(tmpStoich)) | as.numeric(tmpStoich)!=as.integer(tmpStoich)) {
    stop(paste("Detected non-integer stoichiometry in ",reactantOrProduct," id ", tmpId,sep=""))
  }
  tmpStoich <- as.integer(tmpStoich)
  #rxtnt <- c(id=tmpId,stoichiometry=tmpStoich)
  reactantOrProduct <- list(id=tmpId,stoichiometry=tmpStoich)
  return(reactantOrProduct)
}

getReactionInfo <- function(reactionXMLNode) {
  if (!any("Id"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    stop("ERROR: Reaction without Id tag")
  }
  if (!any("Type"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    stop("ERROR: Reaction without Type tag")
  }
  # StochKit2R does not allow Type=customized, therefore Rate is required
  if (XML::xmlValue(reactionXMLNode[["Type"]])!="mass-action") {
    stop("ERROR: Reaction Type must equal 'mass-action'")
  }
  if (!any("Rate"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    stop("ERROR: Reaction without Rate tag")
  }
  #must contain at least one of Reactants or Products (usually both)
  #note, this does not catch case where no SpeciesReference tags are present
  if (!any("Reactants"==XML::xmlSApply(reactionXMLNode, XML::xmlName)) & !any("Products"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    stop("ERROR: Reaction must contain Reactants or Products (or both)")
  }
  
  #process Reactants
  Reactants=list()
  if (any("Reactants"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    # if no SpeciesReference tags, will get an empty $text element
    # so, check for SpeciesReference in the names
    if (sum(names(XML::xmlChildren(reactionXMLNode[["Reactants"]]))=="SpeciesReference")>0) {
      reactantsRaw <- XML::xmlChildren(reactionXMLNode[['Reactants']])
      #iterate over children (reactants), appending to list
      for (SR in XML::xmlChildren(reactionXMLNode[["Reactants"]])) {
        Reactants = c(Reactants,list(getReactantOrProduct(SR,"Reactant")))
      }
    }
  }
  #process Products
  Products=list()
  if (any("Products"==XML::xmlSApply(reactionXMLNode, XML::xmlName))) {
    if (sum(names(XML::xmlChildren(reactionXMLNode[["Products"]]))=="SpeciesReference")>0) {
      productsRaw <- XML::xmlChildren(reactionXMLNode[['Products']])
      #iterate over children (reactants), appending to list
      for (SR in XML::xmlChildren(reactionXMLNode[["Products"]])) {
        Products = c(Products,list(getReactantOrProduct(SR,"Product")))
      }
    }
  }
  
  # catch case where Reactants and Products tags both exist but are empty
  if (length(Products)==0 & length(Reactants)==0) {
    stop("Reactants and Products contain no SpeciesReference entries. At least one reactant or product must be specified")
  }
  
  rxns <- list(Id=XML::xmlValue(reactionXMLNode[["Id"]]),Rate=XML::xmlValue(reactionXMLNode[["Rate"]]),Reactants=Reactants,Products=Products)
  if (sub("[[:space:]]+$", "", rxns["Id"])!=rxns["Id"]) {
    warning(paste("White space detected in Reaction Id ",rxns["Id"],", removing spaces"))
    rxns["Id"]=sub("[[:space:]]+$", "", rxns["Id"])
  }
  return(rxns)
}

