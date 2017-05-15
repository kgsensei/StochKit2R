buildStochKit2Rmodel <- function(modelFile) {
#  doc <- xmlInternalTreeParse("/Users/kevinsanft/Desktop/StochKit2.0.10/models/examples/michaelis_menten_andreas_full.xml")
#  doc <- xmlInternalTreeParse("/Users/kevinsanft/Desktop/testmodel.xml")
  
  doc <- XML::xmlInternalTreeParse(modelFile)
  
  #how many children? #should be 1: model
  if (XML::xmlSize(doc)!=1) {
    stop("xml model file should have one top level tag: Model")
  } 
  model = XML::xmlRoot(doc)
  if (XML::xmlName(model)!="Model") {
    stop("xml model file should have one top level tag: Model")
  }
  # StochKit2R cannot handle Events
  if (length(XML::getNodeSet(doc,"//EventsList"))>0) {
    stop("Events not supported in StochKit2R")
  }
  
  #how many children?
  #xmlSize(model) #6 for dimer decay model: Description, NumberOfReactions, NumberOfSpecies, ParametersList, ReactionsList, SpeciesList
  
  # the children of model contain the important model information
  modelChildren <- XML::xmlChildren(model)
  
  StochKitNumberOfReactions <- as.numeric(XML::xmlValue(modelChildren[["NumberOfReactions"]]))
  if (is.na(StochKitNumberOfReactions)) {
    stop("NumberOfReactions tag missing or value invalid")
  }
  StochKitNumberOfSpecies <- as.numeric(XML::xmlValue(modelChildren[["NumberOfSpecies"]]))
  if (is.na(StochKitNumberOfSpecies)) {
    stop("NumberOfSpecies tag missing or value invalid")
  }
  
  # create StochKit2RParameterList. Empty list of no Parameters
  StochKit2RParameterList=list()
  #first we need to get the Parameters from ParametersList
  StochKitParametersList <- modelChildren[["ParametersList"]]
  # will be NULL of no ParametersList
  if (!is.null(StochKitParametersList)) {
    #all children of ParametersList should have name=="Parameter"
    if (!all("Parameter"==XML::xmlSApply(StochKitParametersList, XML::xmlName))) {
      stop("Unexpected tag in ParametersList. Expected only Parameter tags")
    }
    
    #how many parameters?
    #NumberOfParameters <- xmlSize(StochKitParametersList) #4 for dimer decay
    StochKitParameters <- XML::xmlChildren(StochKitParametersList)
    
    StochKit2RParameterList <- lapply(StochKitParameters,getParameterIdAndExpression)
  }
  
  # next, let's get SpeciesList
  StochKitSpeciesList <- modelChildren[["SpeciesList"]]
  if (is.null(StochKitSpeciesList)) {
    stop("Missing required SpeciesList tag")
  }
  #all children should have name=="Species"
  if (!all("Species"==XML::xmlSApply(StochKitSpeciesList, XML::xmlName))) {
    stop("Unexpected tag in SpeciesList. Expected only Species tags")
  }
  
  #how many species?
  StochKitSpeciesListSize <- XML::xmlSize(StochKitSpeciesList) #3 for dimer decay
  if (StochKitSpeciesListSize!=StochKitNumberOfSpecies) {
    stop("NumberOfSpecies does not match number of species in SpeciesList")
  }
  StochKitSpecies <- XML::xmlChildren(StochKitSpeciesList)
    
  StochKit2RSpeciesList <- lapply(StochKitSpecies,getSpeciesIdAndInitialPopulation)
  
  #ReactionsList
  StochKitReactionsList <- modelChildren[["ReactionsList"]]
  StochKitReactions <- XML::xmlChildren(StochKitReactionsList)
  StochKitReactionsListSize <- XML::xmlSize(StochKitReactions)
  if (StochKitReactionsListSize!=StochKitNumberOfReactions) {
    stop("ERROR: NumberOfReactions does not match number of reactions in ReactionList")
  }
  
  StochKit2RReactionList <- lapply(StochKitReactions,getReactionInfo)
  
  return(list(StochKit2RParameterList,StochKit2RSpeciesList,StochKit2RReactionList))
}