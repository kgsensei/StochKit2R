#include "solver_helper_functions.h"

int assignRealizations(int totalRealizations, int totalProcesses, int id) {
    //assign realizations to threads
    //if totalRealizations%totalProcesses==0, easy
    //but if not, assign the remaining N realizations to the first N ids
    int assignedRealizations=totalRealizations/totalProcesses;
    //assign up to one more realization
    if (id < totalRealizations%totalProcesses) {
        assignedRealizations+=1;
    }
    return assignedRealizations;
}

std::vector<std::string> getSpeciesList(Rcpp::List rSpeciesList) {
    std::vector<std::string> modelSpeciesList;
    for (int i=0; i!=rSpeciesList.size(); ++i) {
        Rcpp::CharacterVector thisSpec = rSpeciesList[i];
        std::string sId = Rcpp::as<std::string>(thisSpec[0]);
        modelSpeciesList.push_back(sId);
    }
    return modelSpeciesList;
}

std::string chooseMethod(int numberOfReactions, int numberOfRealizations, int processes) {
    if (numberOfReactions>constantOverODMCutoff) {
        //large M, use constant-complexity method
        return "constant";
    } else {
        if (numberOfRealizations>=realizationPerThreadCutoff*processes && numberOfReactions>1) {
            //if we're here, we prefer ODM over Direct
            return "odm";
        } else {
            //if we're here, we prefer Direct over ODM
            //see if constant would be better
            if ( (numberOfReactions/constantOverDirectCutoff+numberOfReactions*numberOfRealizations/constantOverODMCutoff) > numberOfRealizations ) {
                //use constant
                return "constant";
            }
            else {
                //use direct
                return "direct";
            }
        }
    }
}
