// [[Rcpp::depends(BH)]]

#include <iostream>
#include <string>
#include "boost/numeric/ublas/vector_sparse.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix_sparse.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "CustomPropensitySet.h"
#include "SSA_Direct.h"
#include "TauLeapingExplicitAdaptive.h"
#include "StandardDriverUtilities.h"
#include <boost/random/mersenne_twister.hpp>
#include "MassActionModelMatrixStoichiometry.h"
#include "solver_helper_functions.h"

#include <Rcpp.h>

#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace Rcpp;

//'@title C++ Interface to Explicit Adaptive Tau-Leaping simulation
//'
//'@description
//'\code{ssa} Called by StochKit2R tauLeaping function, do not call this C++ interface directly
//'
//'@param StochKit2Rmodel R list (Rcpp List built from buildStochKit2Rmodel output)
//'@param outputDirNameString Character string with path to output directory. Should end in path separator.
//'@param time Simulation time of each realization
//'@param realizations Number of realizations
//'@param intervals Number of output intervals. Default 0 outputs at end time only. 1=keep data at start and end time, 2=keep data at start, middle, and end times, etc. Note data is stored at (intervals+1) equally spaced time points.
//'@param keepStats Keep means and variances data
//'@param keepTrajectories Keep trajectory data
//'@param keepHistograms Keep histogram data
//'@param bins Number of histogram bins
//'@param seed Seed for random number generator
//'@param p Override default and specify the number of processes (threads) to use. By default (=0), the number of processes will be determined automatically
//'@param epsilon Set the tolerance (applicable to tauLeaping only), default is 0.03. Valid values: must be greater than 0.0 and less than 1.0
//'@param threshold Set the threshold (minimum number of reactions per leap before switching to ssa) for tauLeaping
//'@return NULL
//'@keywords internal
// [[Rcpp::export]]
void tauLeapingStochKit2RInterface(Rcpp::List StochKit2Rmodel, std::string outputDirNameString, double time, int realizations, int intervals, bool keepStats, bool keepTrajectories, bool keepHistograms, int bins, unsigned int seed, int p, double epsilon, int threshold) {
    //assumes outputDirNameString is a valid path to the output directory name that does not end in path separator

    //create StochKit2R mass action model object
    //first, pull out pieces from list object
    Rcpp::List rParameterList=StochKit2Rmodel[0];
    Rcpp::List rSpeciesList=StochKit2Rmodel[1];
    Rcpp::List rReactionList=StochKit2Rmodel[2];
    
    STOCHKIT::MassActionModelMatrixStoichiometry<STOCHKIT::StandardDriverTypes::populationType, STOCHKIT::StandardDriverTypes::matrixStoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,rReactionList,rSpeciesList);
    
    //create the main output object
    std::vector<STOCHKIT::StandardDriverTypes::outputType> output(1);
    //set output options
    output[0].setOutputTimes(STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::createUniformOutputTimes(0.0,time,intervals));
    output[0].setKeepStats(keepStats);
    output[0].setKeepTrajectories(keepTrajectories);
    output[0].setKeepHistograms(keepHistograms);
    output[0].setHistogramBins(bins);
    
    //use StochKit2 seeding (and parallel RNG) strategy
    std::vector<unsigned int> seeds;
    boost::mt19937 seedGenerator;
    seedGenerator.seed(seed);
    seeds.push_back(seedGenerator());//thread 0
    
    
    Rcout << "StochKit2R MESSAGE: running tau-leaping...\n";

#if defined(_OPENMP)
    int n;//number of threads
    if (p!=0) {
        omp_set_num_threads(p); //force use of p threads per user's request
    }
#pragma omp parallel
    {
        int ID = omp_get_thread_num();
        
#pragma omp single
        {
            n=omp_get_num_threads();
            output.resize(n);
            //seed other threads
            for (int i=1; i<n; ++i) {
                seeds.push_back(seedGenerator());
            }
        }
//ensure all RNGs are seeded before proceeding
#pragma omp barrier
        int localRealizations=assignRealizations(realizations, n, ID);

        if (ID!=0) {
            //set output options
            output[ID].setOutputTimes(STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::createUniformOutputTimes(0.0,time,intervals));
            output[ID].setKeepStats(keepStats);
            output[ID].setKeepTrajectories(keepTrajectories);
            output[ID].setKeepHistograms(keepHistograms);
            output[ID].setHistogramBins(bins);
        }
        
        STOCHKIT::TauLeapingExplicitAdaptive<STOCHKIT::StandardDriverTypes::populationType, STOCHKIT::StandardDriverTypes::matrixStoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> tauLeaping(model.writeInitialPopulation(),model.writeMatrixStoichiometry(),model.writePropensities(),model.writeDependencyGraphMatrixStoichiometry(),seeds[ID]);
        tauLeaping.setEpsilon(epsilon);
        tauLeaping.setThreshold(threshold);
        tauLeaping.simulate<STOCHKIT::StandardDriverTypes::outputType>(localRealizations, 0.0, time, output[ID]);

        //merge output
        //in the future this can be done in parallel and use a little less data
#pragma omp master
        {
            for (int i=1; i<n; ++i) {
                output[0].merge(output[i]);
            }

        }
    }
#else
    STOCHKIT::TauLeapingExplicitAdaptive<STOCHKIT::StandardDriverTypes::populationType, STOCHKIT::StandardDriverTypes::matrixStoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> tauLeaping(model.writeInitialPopulation(),model.writeMatrixStoichiometry(),model.writePropensities(),model.writeDependencyGraphMatrixStoichiometry(),seeds[0]);
    tauLeaping.setEpsilon(epsilon);
    tauLeaping.setThreshold(threshold);
     tauLeaping.simulate<STOCHKIT::StandardDriverTypes::outputType>(realizations, 0.0, time, output[0]);

#endif
    
    //write output
    //get species labels...
    std::vector<std::string> modelSpeciesList=getSpeciesList(rSpeciesList);
    
    if (keepStats) {
        Rcout << "creating statistics output files...\n";
        STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/means.txt",modelSpeciesList);
        output[0].stats.writeMeansToFile(outputDirNameString+"/stats/means.txt",true,true);
        STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/variances.txt",modelSpeciesList);
        output[0].stats.writeVariancesToFile(outputDirNameString+"/stats/variances.txt",true,true);
    }
    if (keepTrajectories) {
        Rcout << "creating trajectories output files...\n";
        std::size_t trajectoryNumber;
        std::string trajectoryNumberString;
        for (int i=0; i!=realizations; ++i) {
            trajectoryNumber=i;
            trajectoryNumberString=STOCHKIT::StandardDriverUtilities::size_t2string(trajectoryNumber);

            //
            STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",modelSpeciesList);
            output[0].trajectories.writeDataToFile(i,outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",true,true);
        }
    }
    if (keepHistograms) {
        Rcout << "creating histogram output files...\n";
        output[0].histograms.writeHistogramsToFile(outputDirNameString+"/histograms/hist",".dat",modelSpeciesList);
    }
    Rcout << "done!\n";
}
