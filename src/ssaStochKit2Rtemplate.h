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
Rcpp::List ssaStochKit2Rtemplate(Rcpp::List& StochKit2Rmodel, double time, int realizations, int intervals, bool keepStats, bool keepTrajectories, bool keepHistograms, int bins, std::string outputDirNameString, unsigned int seed, int p) {
  //assumes outputDirNameString is a valid path to the output directory name that does not end in path separator

  //create StochKit2R mass action model object
  //first, pull out pieces from list object
  Rcpp::List rParameterList=StochKit2Rmodel[0];
  Rcpp::List rSpeciesList=StochKit2Rmodel[1];
  Rcpp::List rReactionList=StochKit2Rmodel[2];
  Rcpp::List rCustomPropensityList=StochKit2Rmodel[3];
//    int NumberOfSpecies=rSpeciesList.size();
//    int NumberOfReactions=rReactionList.size();
    
  STOCHKIT::MassActionModel<STOCHKIT::StandardDriverTypes::populationType, _stoichiometryType,
    STOCHKIT::StandardDriverTypes::propensitiesType, 
	STOCHKIT::StandardDriverTypes::graphType> model(rParameterList,rReactionList,rSpeciesList,rCustomPropensityList);
    
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
	
	int n=1; //number of threads, may be set >1 later
	int defaultN=1;//default number of threads, so we can reset later
	
    
#if defined(_OPENMP)
#pragma omp parallel
	{
#pragma omp master
		{
			defaultN = omp_get_num_threads();//so we can reset later
		}
	}

  if (p!=0) {
    omp_set_num_threads(p); //force use of p threads per user's request
  }
#pragma omp parallel
  {
    int ID = omp_get_thread_num();
        
#pragma omp single
    {
      n=omp_get_num_threads();
      Rcpp::Rcout << "...running "<<n<<" threads...\n";
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
			Rcout << "...merging output...\n";
            for (int i=1; i<n; ++i) {
                output[0].merge(output[i]);
            }
        }
    }
	//reset number of threads to default
	omp_set_num_threads(defaultN);

#else
    
  // we are in serial mode (no OpenMP support)
  Rcpp::Rcout << "...running in serial mode (1 thread)...\n";
  _solverType ssa(model.writeInitialPopulation(),model.writeStoichiometry(),model.writePropensities(),model.writeDependencyGraph(),seeds[0]);
  ssa.template simulate<STOCHKIT::StandardDriverTypes::outputType>(realizations, 0.0, time, output[0]);

#endif
    
  //simulation is complete, all output is stored in output[0]
  //write output to files
  //get species labels...
  std::vector<std::string> modelSpeciesList=getSpeciesList(rSpeciesList);

  int nr, nc;
  // name the columns accordingly
  Rcpp::CharacterVector col_names;
  col_names.push_back("time");

  for (int i = 0; i < modelSpeciesList.size(); i++)
    col_names.push_back(modelSpeciesList[i]);

  Rcpp::DataFrame means = R_NilValue;
  Rcpp::DataFrame vars = R_NilValue;
  Rcpp::List Stats = R_NilValue;
  Rcpp::List trajs = R_NilValue;
  Rcpp::List hist = R_NilValue;
	
  if (keepStats) {
    if(outputDirNameString != "") {
      Rcpp::Rcout << "creating statistics output files...\n";
      STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/means.txt",modelSpeciesList);
      output[0].stats.writeMeansToFile(outputDirNameString+"/stats/means.txt",true,true);
      STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/stats/variances.txt",modelSpeciesList);
      output[0].stats.writeVariancesToFile(outputDirNameString+"/stats/variances.txt",true,true);
    }
    // return data
    // each 2D double array will get converted to a dataframe
    // retrieve stats in 2D arrays
    std::vector< std::vector<double> > meansBuffer = output[0].stats.getMeans();
    std::vector< std::vector<double> > variancesBuffer = output[0].stats.getVariances();

    //////////// MEANS ///////////////
    nr = meansBuffer.size();
    nc = meansBuffer[0].size();

    Rcpp::NumericMatrix m( nr, nc );
    for (int i=0; i<nr; i++) {
      std::vector<double>& result_i = meansBuffer[i];

      if (result_i.size() != nc) //the length of data in a row is inconsistent
        Rcpp::stop("Fatal error encountered, terminating StochKit2R"); //cannot convert to data frame.

      for (int j=0; j<nc; j++)
        m(i, j) = result_i[j];
    }

    means = Rcpp::DataFrame(m);
    means.attr("names")= col_names;
    ////////////////////////////

    //////////// VARIANCES ///////////////
    nr = variancesBuffer.size();
    nc = variancesBuffer[0].size();

    Rcpp::NumericMatrix v( nr, nc );
    for (int i=0; i<nr; i++) {
      std::vector<double>& result_i = variancesBuffer[i];

      if (result_i.size() != nc) //the length of data in a row is inconsistent
        Rcpp::stop("Fatal error encountered, terminating StochKit2R"); //cannot convert to data frame.

      for (int j=0; j<nc; j++)
        v(i, j) = result_i[j];
    }

    vars = Rcpp::DataFrame(v);
    vars.attr("names")= col_names;
    ///////////////////////////
    }


  if (keepTrajectories) {
    if(outputDirNameString!="") {
      Rcpp::Rcout << "creating trajectories output files...\n";
      std::size_t trajectoryNumber;
      std::string trajectoryNumberString;
        for (int i=0; i!=realizations; ++i) {
          trajectoryNumber=i;
          trajectoryNumberString=STOCHKIT::StandardDriverUtilities::size_t2string(trajectoryNumber);
      
          STOCHKIT::IntervalOutput<STOCHKIT::StandardDriverTypes::populationType>::writeLabelsToFile(outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",modelSpeciesList);
          output[0].trajectories.writeDataToFile(i,outputDirNameString+"/trajectories/trajectory"+trajectoryNumberString+".txt",true,true);
      }
    }
    // retrieve trajs in 3D array
    std::vector< std::vector< std::vector<double> > > trajBuffer;

    for (int i=0; i!=realizations; ++i)
      trajBuffer.push_back( output[0].trajectories.getData(i) );

    //////////// TRAJECTORIES ///////////////
    trajs = Rcpp::List(trajBuffer.size());

    for ( int ii = 0; ii < trajBuffer.size(); ii++ ) {
      nr = trajBuffer[ii].size();
      nc = trajBuffer[ii][0].size();

      Rcpp::NumericMatrix st( nr, nc );
      for (int i=0; i<nr; i++) {
        std::vector<double>& result_i = trajBuffer[ii][i];

        if (result_i.size() != nc) //the length of data in a row is inconsistent
          Rcpp::stop("Fatal error encountered, terminating StochKit2R");//cannot convert to data frame

        for (int j=0; j<nc; j++)
          st(i, j) = result_i[j];
      }

      Rcpp::DataFrame singTraj = Rcpp::DataFrame(st);
      singTraj.attr("names")= col_names;

      trajs[ii] = singTraj;
    }
  }

  if (keepHistograms) {
  	if(outputDirNameString!=""){  
      Rcpp::Rcout << "creating histogram output files...\n";
  	  output[0].histograms.writeHistogramsToFile(outputDirNameString+"/histograms/hist",".dat",modelSpeciesList);
  	} 
	  // create hist return list
	  // hist has one element per species
	  hist = Rcpp::List(output[0].histograms.numberOfSpecies());
	  // name elements species names
	  hist.attr("names")= modelSpeciesList;
	  
	  // iterate over the species
	  for (int species_index = 0; species_index<modelSpeciesList.size(); species_index++) {

		  std::vector<std::vector<std::string> > species_data;
		  //iterate over the output intervals
		  for (int output_index = 0; output_index<=intervals; output_index++) {
			  std::vector<std::string> this_interval = output[0].histograms(output_index,species_index).fileDataAsString(output[0].histograms.getOutputTimes(), modelSpeciesList[species_index]);
			  species_data.push_back(this_interval);
		  }
		  
		  hist[species_index] = Rcpp::wrap(species_data);
	  }
  }

  //merge means and vars into "stats"
  Rcpp::List stats = Rcpp::List::create(
    Rcpp::Named("means") = means,
    Rcpp::Named("vars") = vars);


  return Rcpp::List::create(
              Rcpp::Named("stats") = stats,
              Rcpp::Named("trajs") = trajs,
							Rcpp::Named("hist") = hist);
}
#endif
