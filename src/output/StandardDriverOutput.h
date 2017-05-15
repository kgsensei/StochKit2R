/******************************************************************************
 *  FILE:    StandardDriverOutput.h
 */

#ifndef _STANDARD_DRIVER_OUTPUT_H_
#define _STANDARD_DRIVER_OUTPUT_H_

#include <iostream>
#include <vector>
#include "IntervalOutput.h"
#include "StatsOutput.h"
#include "HistogramOutput.h"
//#include "StandardDriverUtilities.h"

namespace STOCHKIT
{
 template<typename _populationVectorType>
  class StandardDriverOutput
 {
public:

 public:
  IntervalOutput<_populationVectorType> trajectories;
  StatsOutput<_populationVectorType> stats;
  HistogramOutput<_populationVectorType> histograms;

 protected:
  bool keepStats;
  bool keepTrajectories;
  bool keepHistograms;

 public:

  StandardDriverOutput(bool keepStats=true,bool keepTrajectories=true, bool keepHistograms=true):trajectories(),stats(),keepStats(keepStats),keepTrajectories(keepTrajectories), keepHistograms(keepHistograms)
  {}
  
  void setKeepStats(bool keep) {
    keepStats=keep;
  }
  void setKeepTrajectories(bool keep) {
    keepTrajectories=keep;
  }
  void setKeepHistograms(bool keep) {
    keepHistograms=keep;
  }

  void setOutputTimes(std::vector<double> outputTimes){
    trajectories.setOutputTimes(outputTimes);
    stats.setOutputTimes(outputTimes);
    histograms.setOutputTimes(outputTimes);
  }

  std::vector<double> getOutputTimes() {
    return trajectories.getOutputTimes();
  }

  bool initialize(std::size_t realizations, double startTime, double endTime, _populationVectorType& samplePopulationVector) {

    bool statsInitialized=true;
    bool trajectoriesInitialized=true;
    bool histogramsInitialized=true;
    //initialize stats
    if (keepStats) {
      statsInitialized=stats.initialize(realizations,startTime, endTime, samplePopulationVector);
    }
    //initialize trajectory data
    if (keepTrajectories) {
	  trajectoriesInitialized=trajectories.initialize(realizations, startTime, endTime, samplePopulationVector);
    }
    //initialize histogram data
    if (keepHistograms) {
      histogramsInitialized=histograms.initialize(realizations, startTime, endTime, samplePopulationVector);
    }

    if (statsInitialized && trajectoriesInitialized && histogramsInitialized) {
      return true;
    }
    else {
      //should provide an error message
      Rcpp::Rcout << "StochKit ERROR (StandardDriverOutput::initialize) initialization of output object failed" << std::endl;
      return false;
    }
  }
  
  void setSpeciesSubset(std::vector<std::size_t> speciesIndices) {
    trajectories.setSpeciesSubset(speciesIndices);
    stats.setSpeciesSubset(speciesIndices);
    histograms.setSpeciesSubset(speciesIndices);
  }
  
  void setHistogramBins(std::size_t bins) {
    histograms.setNumberOfBins(bins);
  }
  
  void record(std::size_t realization, std::size_t interval,_populationVectorType population) {
    if (keepTrajectories) {
      trajectories.record(realization,interval,population);
    }
    if (keepStats) {
      stats.record(realization,interval,population);
    }
    if (keepHistograms) {
      histograms.record(realization,interval,population);
    }
  }
  
  void writeMeansToFile(std::string filename, bool printTime=true, bool append=false) {
    stats.writeMeansToFile(filename,printTime,append);
  }
  
  void writeVariancesToFile(std::string filename, bool printTime=true, bool append=false) {
    stats.writeVariancesToFile(filename,printTime,append);
  }
  
  void writeStandardDeviationsToFile(std::string filename, bool printTime=true, bool append=false) {
    stats.writeStandardDeviationsToFile(filename,printTime,append);
  }
  
  void writeTrajectoryToFile(std::size_t realization, std::string filename, bool printTime=true, bool append=false) {
    if (keepTrajectories) {
      trajectories.writeDataToFile(realization,filename,printTime,append);
    }
    else {
      Rcpp::Rcout << "WARNING (StandardDriverOutput::writeTrajectoryToFile): keepTrajectories=false, can't write trajectory data to file" << std::endl;
    }
  }

  void merge(StandardDriverOutput<_populationVectorType> other) {
    stats.merge(other.stats);
    trajectories.merge(other.trajectories);
  }

  void merge(std::vector<StandardDriverOutput<_populationVectorType> > others) {
    for (std::size_t i=0; i!=others.size(); ++i) {
      this->merge(others[i]);
    }
  }

 };
}

#endif
  
