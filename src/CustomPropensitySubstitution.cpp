#include "Parameter.h"
#include <Rcpp.h>
using namespace Rcpp;

class TempSpecies{
public:
	std::string Id; // Id of the species
	std::string InitialPopulation;  // initial population, only valid in SpeciesList
	std::vector<int> AffectReactions; // only valid in SpeciesList, record what reactions does one Species affect
};

std::string customPropensitySubstitution(std::string equation, STOCHKIT::ListOfParameters &ParametersList, std::vector<TempSpecies> &SpeciesList)
{
	// locate a parameter name in equation
	unsigned int begin = 0, end = 0;
	std::string substitutedEquation = equation;
	
	while( begin < substitutedEquation.length() ){
		if( (isalpha(substitutedEquation.at(begin)) || substitutedEquation.at(begin) == '_') && ( (begin==0) || ((substitutedEquation.at(begin) != 'e') && (substitutedEquation.at(begin) != 'E')) || !(isalnum(substitutedEquation.at(begin-1)) || substitutedEquation.at(begin-1) == '_') )){
			end = begin+1;
			while( (end<substitutedEquation.length()) && (isalnum(substitutedEquation.at(end)) || substitutedEquation.at(end) == '_') )
				++end;
			
			std::string parameterName = substitutedEquation.substr(begin,end-begin);
			
			// search for the parameter in ParametersList
			unsigned int i = 0;
			while( (i < ParametersList.size()) && (parameterName.compare(ParametersList[i].Id)!=0) )
				++i;
			
			if( i != ParametersList.size() ){
				// substitute parameter with its value
				std::ostringstream parameterValue;
				parameterValue <<  ParametersList[i].Value;
				substitutedEquation.replace(begin, end-begin, parameterValue.str());
				
				begin = begin + parameterValue.str().size();
			}
			else{
				// search for species in SpeciessList
				unsigned int j = 0;
				while( (j < SpeciesList.size()) && (parameterName.compare(SpeciesList[j].Id)!=0) )
					++j;
				
				if( j != SpeciesList.size() ){
					// substitute species name with its population variable
					std::ostringstream speciesReference;
					speciesReference << "x[" << j <<"]";
					substitutedEquation.replace(begin, end-begin, speciesReference.str());
					
					begin = begin + speciesReference.str().size();
				}
				else{
					//std::cout << "StochKit NOTE (CustomPropensitySubstitution): function \"" << parameterName << "\" written into custom propensity function (compilation will fail if this is not a valid C++ function).\n";
					begin = begin + parameterName.size();
				}
			}
		}
		else if (substitutedEquation.at(begin) == '/'){
			substitutedEquation.replace(begin, 1, "/(double)");
			begin += 9;
		}
		else{
			++begin;
		}
	}
	
	return substitutedEquation;
}

// [[Rcpp::export]]
StringVector customPropensitySubstitution(Rcpp::StringVector originalString, List originalParametersList, List originalSpeciesList) {

	//return originalParametersList[0];
	//create ListOfParameters object from originalParametersList
	STOCHKIT::Parameter *cur_ptr = NULL;
	STOCHKIT::ListOfParameters listOfParameters;
	for (int i = 0; i!=originalParametersList.size(); i++) {
		listOfParameters.push_back(STOCHKIT::Parameter());
		cur_ptr = &listOfParameters.back(); // reference to current parameter
		
		StringVector rParam = originalParametersList[i];
		
		
		cur_ptr->Id = rParam[0];
		cur_ptr->Expression = rParam[1];
		cur_ptr->Type = 0;
		cur_ptr->CalculateFlag = -1;

		//Rcpp::Rcout << cur_ptr->Id << ", "<<cur_ptr->Expression<<"\n";
	}
	
	if(!listOfParameters.linkParameters()){
		Rcpp::Rcout << "StochKit ERROR (customPropensitySubstitution): error linking parameters.\n";
		Rcpp::stop("Fatal error encountered, terminating StochKit2R");
	}

	if(!listOfParameters.calculateParameters()){
		Rcpp::Rcout << "StochKit ERROR (customPropensitySubstitution): error calculating parameters.\n";
		Rcpp::stop("Fatal error encountered, terminating StochKit2R");
	}

//	if(!linkSpeciesAndReactions()){
//		exit(1);
//	}
	
	//need SpeciesList duplicated code :( named TempSpecies
	std::vector<TempSpecies> SpeciesList;
	TempSpecies *cur_ptr2 = NULL;

	for (int i = 0; i!=originalSpeciesList.size(); i++) {

		SpeciesList.push_back(TempSpecies());
		cur_ptr2 = &SpeciesList.back();
		StringVector rSpecies = originalSpeciesList[i];
		
		
		cur_ptr2->Id = rSpecies[0];
		cur_ptr2->InitialPopulation = rSpecies[1];

	}
	
	//create output vector
	StringVector output(originalString.size());
	//iterate over input strings, convert, put in output
	for (int i =0; i<originalString.size(); i++) {
		output[i] = NA_STRING;//customPropensitySubstitution(str, listOfParameters, SpeciesList);
		//if string is not NA, convert it
		if (!StringVector::is_na(originalString[i]) ) {
			std::string str = as<std::string>(originalString[i]);
			output[i]=customPropensitySubstitution(str, listOfParameters,SpeciesList);
		}
	}

	
	return output;//result;//originalString;//"test123";
}


