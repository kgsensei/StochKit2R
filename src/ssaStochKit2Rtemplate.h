#ifndef SSASTOCHKIT2RTEMPLATE_H
#define SSASTOCHKIT2RTEMPLATE_H

// [[Rcpp::depends(BH)]]

#include <iostream>
#include <string>
#include "boost/numeric/ublas/vector_sparse.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix_sparse.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include "CustomPropensitySet.h"
#include "SSA_Direct.h"
#include "SSA_ODM.h"
#include "SSA_ConstantTime.h"
#include "StandardDriverUtilities.h"
#include <boost/random/mersenne_twister.hpp>
#include "MassActionModel.h"
#include "solver_helper_functions.h"

#include <Rcpp.h>

#if defined(_OPENMP)
#include "omp.h"
#endif

using namespace Rcpp;

template<typename _stoichiometryType,
typename _solverType>
void ssaStochKit2Rtemplate(Rcpp::List& StochKit2Rmodel, std::string outputDirNameString, double time, int realizations, int intervals, bool keepStats, bool keepTrajectories, bool keepHistograms, int bins, unsigned int seed, int p) {
    //assumes outputDirNameString is a valid path to the output directory name that does not end in path separator

    //create StochKit2R mass action model object
    //first, pull out pieces from list object
    Rcpp::List rParameterList=StochKit2Rmodel[0];
    Rcpp::List rSpeciesList=StochKit2Rmodel[1];
    Rcpp::List rReactionList=StochKit2Rmodel[2];
    
//    int NumberOfSpecies=rSpeciesList.size();
//    int NumberOfReactions=rReactionList.size();
    
    STOCHKIT::MassActionModel<STOCHKIT::StandardDriverTypes::populationType, _stoichiometryType, STOCHKIT::StandardDriverTypes::propensitiesType, STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,rReactionList,rSpeciesList);
    
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
    
    
#if defined(_OPENMP)
    int n=1; //number of threads, may be set >1 later

    if (p!=0) {
        omp_set_num_threads(p); //force use of p threads per user's request
    }
#pragma omp parallel
    {
        int ID = omp_get_thread_num();
        
#pragma omp single
        {
            n=omp_get_num_threads();
            Rcout << "...running "<<n<<" threads...\n";
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
        
        _solverType ssa(model.writeInitialPopulation(),model.writeStoichiometry(),model.writePropensities(),model.writeDependencyGraph(),seeds[ID]);
        ssa.template simulate<STOCHKIT::StandardDriverTypes::outputType>(localRealizations, 0.0, time, output[ID]);
        
        //merge output
        //in the future this can be done in parallel and use a little less data
#pragma omp barrier
#pragma omp master
        {
            for (int i=1; i<n; ++i) {
                output[0].merge(output[i]);
            }
        }
    }
#else
    
    // we are in serial mode (no OpenMP support)
    Rcout << "...running in serial mode (1 thread)...\n";
    _solverType ssa(model.writeInitialPopulation(),model.writeStoichiometry(),model.writePropensities(),model.writeDependencyGraph(),seeds[0]);
    ssa.template simulate<STOCHKIT::StandardDriverTypes::outputType>(realizations, 0.0, time, output[0]);


#endif
    
    //simulation is complete, all output is stored in output[0]
    //write output to files
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

}
#endif
