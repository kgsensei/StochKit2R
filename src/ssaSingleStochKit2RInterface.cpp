// [[Rcpp::depends(BH)]]

#include <iostream>
#include <fstream>
#include <string>
#include "boost/numeric/ublas/vector_sparse.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix_sparse.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "CustomPropensitySet.h"
#include "SSA_Direct.h"
#include "StandardDriverUtilities.h"
#include <boost/random/mersenne_twister.hpp>
#include "MassActionModel.h"
#include "solver_helper_functions.h"

#include <Rcpp.h>

using namespace Rcpp;

//'@title C++ Interface to Gillespie Stochastic Simulation Algorithm single trajectory
//'
//'@description
//'\code{ssa} Called by StochKit2R ssaSingle function, do not call this C++ interface directly
//'
//'@param StochKit2Rmodel R list (Rcpp List built from buildStochKit2Rmodel output)
//'@param outputFileNameString Character string with path to output file.
//'@param startTime Simulation start time
//'@param endTime Simulation end time
//'@param seed Seed for random number generator
//'@return NULL
//'@keywords internal
// [[Rcpp::export]]
void ssaSingleStochKit2RInterface(Rcpp::List StochKit2Rmodel, std::string outputFileNameString, double startTime, double endTime, unsigned int seed) {

    //create StochKit2R mass action model object
    //first, pull out pieces from list object
    Rcpp::List rParameterList=StochKit2Rmodel[0];
    Rcpp::List rSpeciesList=StochKit2Rmodel[1];
    Rcpp::List rReactionList=StochKit2Rmodel[2];
    
    //get species labels...
    std::vector<std::string> modelSpeciesList=getSpeciesList(rSpeciesList);
    //try to write labels to file
    STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputFileNameString,modelSpeciesList);
    
    STOCHKIT::MassActionModel<STOCHKIT::StandardDriverTypes::populationType, STOCHKIT::StandardDriverTypes::denseStoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,rReactionList,rSpeciesList);
    
    //create the output object
    std::vector<std::pair<double, STOCHKIT::StandardDriverTypes::populationType> > output;
    
    //use StochKit2 seeding (and parallel RNG) strategy
    std::vector<unsigned int> seeds;
    boost::mt19937 seedGenerator;
    seedGenerator.seed(seed);
    seeds.push_back(seedGenerator());//thread 0
    
    STOCHKIT::SSA_Direct<STOCHKIT::StandardDriverTypes::populationType, STOCHKIT::StandardDriverTypes::denseStoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> ssa(model.writeInitialPopulation(),model.writeStoichiometry(),model.writePropensities(),model.writeDependencyGraph(),seeds[0]);
     ssa.simulateSingle(startTime, endTime, output);
    
    //write output
    std::ofstream outfile;
    
    //open for appending
    outfile.open(outputFileNameString.c_str(),std::ios::out | std::ios::app);
    
    if (!outfile) {
        Rcpp::Rcout << "StochKit ERROR (ssaSingleStochKit2RInterface): Unable to open output file.\n";
        Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }
    
    try {
        for (std::size_t step=0; step!=output.size(); ++step) {
            //write time
            outfile << output[step].first << "\t";
            
            for (size_t index=0; index!=output[step].second.size(); ++index) {
                outfile << std::setprecision(8) << output[step].second[index] << "\t";
            }
            outfile << "\n";
        }
        outfile.close();
    }
    catch (...) {
        Rcpp::Rcout << "StochKit ERROR (IntervalOutput::writeDataToFile): error writing data to output file.\n";
        Rcpp::stop("Fatal error encountered, terminating StochKit2R");
    }

    
//    if (keepStats) {
//        STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/means.txt",modelSpeciesList);
//        output[0].stats.writeMeansToFile(outputDirNameString+"/stats/means.txt",true,true);
//        STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/variances.txt",modelSpeciesList);
//        output[0].stats.writeVariancesToFile(outputDirNameString+"/stats/variances.txt",true,true);
//    }
//    if (keepTrajectories) {
//        std::size_t trajectoryNumber;
//        std::string trajectoryNumberString;
//        for (int i=0; i!=realizations; ++i) {
//            trajectoryNumber=i;
//            trajectoryNumberString=STOCHKIT::StandardDriverUtilities::size_t2string(trajectoryNumber);
//
//            //
//            STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",modelSpeciesList);
//            output[0].trajectories.writeDataToFile(i,outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",true,true);
//        }
//    }
//    if (keepHistograms) {
//        output[0].histograms.writeHistogramsToFile(outputDirNameString+"/histograms/hist",".dat",modelSpeciesList);
//    }

}
