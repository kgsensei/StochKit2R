createFunctionList <- function(expressionStrings) {
  #take a list of C++ expressions, wrap in C++ functions, compile
  # return a list of C++ functions (or function pointers?)
  
  if (!all(is.na(expressionStrings))) {
    message("NOTE: compiling customized propensity functions, this may take a few moments.")
  }
  
  output = list()
  cppincludes = "#include <boost/numeric/ublas/vector.hpp>"
  functionprefix = "double f"
  functionargs = "(const boost::numeric::ublas::vector<double>& x) {\n\treturn "
  
  functionsuffix = ";\n}"
  for (i in 1:length(expressionStrings)) {
    if (is.na(expressionStrings[i])) {
      output[[i]]=NA
    }
    else {
      cppfun <- paste(functionprefix,i,functionargs,expressionStrings[i],functionsuffix,sep="")#"double h1(boost::numeric::ublas::vector<double>& x1) {return x1(1)*x1(2);}"
      #print(cppfun)
      output[[i]] = cppXPtr(cppfun,depends="BH",includes = cppincludes)
    }
  }
  
  if (!all(is.na(expressionStrings))) {
    message("NOTE: finished compiling customized propensity functions.")
  }
  
  return(output)
}

buildStochKit2Rmodel <- function(modelFile) {
  # modelFile = "~/Desktop/dimer_decay.xml"  
  # modelFile = "~/Desktop/schlogl.xml"
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
  
  # verify that species in ReactionList all appear in SpeciesList
  # build vector of speciesIds
  speciesIds = unlist(lapply(StochKit2RSpeciesList,function(x) return(x["Id"])))
  # check that all reactants and products are in speciesIds
  for (i in 1:length(StochKit2RReactionList)) {
    thisrxn = StochKit2RReactionList[[i]]
    # iterate over reactants
    for (j in thisrxn$Reactants) {
      #print(j$id)
      if (! (j$id %in% speciesIds)) {
        stop(paste("ERROR: species",j$id,"in reaction",thisrxn$Id,"is not in SpeciesList"))
      }
    }
    #iterate over products
    for (j in thisrxn$Products) {
      #print(j$id)
      if (! (j$id %in% speciesIds)) {
        stop(paste("ERROR: species",j$id,"in reaction",thisrxn$Id,"is not in SpeciesList"))
      }
    }    
  }
  
  #get custom PropensityFunction strings
  propFunctionStrings = unlist(lapply(StochKit2RReactionList,function(x) return(x["PropensityFunction"])))
  #create C++ functions out of propensity function strings
  #substitute parameter names for values and species names for 
  expressionStrings = customPropensitySubstitution(propFunctionStrings,StochKit2RParameterList,StochKit2RSpeciesList)
  #convert expressionStrings into list of functions or function pointers
  StochKit2CustomPropensityFunctionList = createFunctionList(expressionStrings)
  #return(list(StochKit2RParameterList,StochKit2RSpeciesList,StochKit2RReactionList))  
  return(list(StochKit2RParameterList,StochKit2RSpeciesList,StochKit2RReactionList,StochKit2CustomPropensityFunctionList))
}