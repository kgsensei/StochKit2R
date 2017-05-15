#ifndef SOLVERHELPERFUNCTIONS_H
#define SOLVERHELPERFUNCTIONS_H

#include <vector>
#include <string>
#include <Rcpp.h>

const int constantOverODMCutoff=2000;
const int constantOverDirectCutoff=1000;
const int realizationPerThreadCutoff=10;
const int denseStoichiometryCutoff=64;//if N < this, use dense stoichiometry vectors

int assignRealizations(int totalRealizations, int totalProcesses, int id);

std::vector<std::string> getSpeciesList(Rcpp::List rSpeciesList);

std::string chooseMethod(int numberOfReactions, int numberOfRealizations, int processes);

#endif
