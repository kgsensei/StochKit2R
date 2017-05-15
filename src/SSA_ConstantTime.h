/*
\brief the amortized constant-complexity (composition-rejection) method of the Stochastic Simulation Algorithm (SSA)
*/

#ifndef _SSA_CONSTANTTIME_H_
#define _SSA_CONSTANTTIME_H_

#include <iostream>
#include <vector>
#include <deque>
#include <list>
#include <utility>
#include <math.h>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "ConstantTimeGroup.h"
#include "ConstantTimeGroupCollection.h"
#include <Rcpp.h>

/*
\file SSA_ConstantTime.h

\brief Gillespie's Stochastic Simulation Algorithm (SSA): constant complexity algorithm

This algorithm is generally faster than the direct method only when the number of reactions is
around 3000-4000 or more.  For a description of the algorithm see
A. Slepoy, A.P. Thompson, and S.J. Plimpton. J Chem Phys 128(20):205101 2008. or
S. Mauch, M. Stalzer "Efficient formulations for exact stochastic simulation of chemical systems"
IEEE/ACM Trans. on Comp. Bio. and Bioinformatics, 30 April 2009.

\param _populationVectorType the population vector type, should be dense
\param _stoichiometryType
\param _propensitiesFunctorType functor takes reaction index and _populationVectorType population
and returns the propensity for that reaction
\param _dependencyGraphType [expand]
*/

namespace STOCHKIT
{
	template<typename _populationVectorType, 
		typename _stoichiometryType,
		typename _propensitiesFunctorType,
		typename _dependencyGraphType>
	class SSA_ConstantTime
	{	
	public:
		typedef _populationVectorType populationVectorType;
		typedef _stoichiometryType stoichiometryType;
		typedef _propensitiesFunctorType propensitiesType;
		typedef _dependencyGraphType dependencyGraphType;

	protected:
		//! the class that implements all random number generator function
		/*! change the RandomGenerator class to swap generators */
		//STOCHKIT::RandomGenerator randomGenerator;
        Random randomGenerator;
        
		//! the initial population
		/*! currentPopulation should be set to initialPopulation at the beginning of each
		realization as in initialize()
		\see initialize()
		*/
		_populationVectorType initialPopulation;
		//! the stoichiometric matrix
		/*!
		should actually be a dense vector of (usually sparse) vectors of dimension NumberOfReactions x NumberOfSpecies
		so that currentPopulation+=stoichiometry[reactionNumber] modifies the population based on
		the stoichiometry of the given reaction number
		*/
		_stoichiometryType stoichiometry;
		//! the propensities functor
		/*! propensities(rxn, pop) returns the propensity of reaction number rxn based on population pop
		after a simulation step, currentPropensities[rxn] should equal propensities(rxn, currentPopulation)
		but since propensities() is a function call, accessing current propensities should be done with
		currentPropensities
		\see currentPropensities
		*/
		_propensitiesFunctorType propensities;
		//! the dependency graph which describes the propensities that are affected by each reaction
		/*!
		should be a dense vector of (usually sparse) vectors of dimension NumberOfReactions x (variable) number of affected reactions
		where dependencyGraph[rxn] returns the variable length vector of reaction indices that are affected by reaction rxn
		e.g. if dependencyGraph[4][0]=2 and dependencyGraph[4][1]=4 and dependencyGraph[4][2]=5 then dependencyGraph[4].size()=3 and
		reaction 4 affects reactions 2, 4, and 5
		\see fireReaction
		*/
		_dependencyGraphType dependencyGraph;

		//! number of species in the system
		std::size_t NumberOfSpecies;
		//! number of reactions in the system
		std::size_t NumberOfReactions;

		//! current time of the simulation, should be incremented at each simulation time step
		double currentTime;
		//! current population of the simulation
		_populationVectorType currentPopulation;
		//! current propensities of the simulation
		std::vector<double> currentPropensities;

		//! index of the last reaction that fired
		/*!
		default and error value is -1
		*/
		int previousReactionIndex;

		//groups of propensities
		ConstantTimeGroupCollection groups;

		//!propensities based on initial population
		std::vector<double> initialPropensities;

		//! groups based on initial conditions
		ConstantTimeGroupCollection initialGroups;


	private:
		//! default constructor not implemented
		SSA_ConstantTime();

	public:

		SSA_ConstantTime(const _populationVectorType& initialPop,
			const _stoichiometryType& stoich,
			const _propensitiesFunctorType& propensitiesFunctor,
			const _dependencyGraphType& depGraph,
			int seed=time(NULL));

		//! destructor
		virtual ~SSA_ConstantTime() {
		}

		/*!
		\brief seed the random number generator
		*/
		void seed(int seed);

		//****these functions were protected and had to be moved to run from single trajectory
		void initialize(double startTime=0.0);
		int selectReaction();
		bool fireReaction(int reactionIndex);
		bool setCurrentTime(double newCurrentTime);
		double selectStepSize();

	protected:
		//! set initial population (not current population), should probably never be used since no 0-argument constructor
		void setInitialPopulation(const _populationVectorType& initialPop) {
			initialPopulation = initialPop;
		}
		//! set stoichiometry matrix, should probably never be used since no 0-argument constructor
		void setStoichiometry(const _stoichiometryType& stoich) {
			stoichiometry=stoich;
		}
		//! set propensities functor, should probably never be used since no 0-argument constructor
		void setPropensities(const _propensitiesFunctorType& propensitiesFunctor) {
			propensities=propensitiesFunctor;
		}
		//! set dependencyGraph, should probably never be used since no 0-argument constructor
		void setDependencyGraph(const _dependencyGraphType& depGraph) {
			dependencyGraph=depGraph;
		}

		/*!
		\brief initialize the state for a new simulation realization, this should be called before each realization
		*/
		

		/*!
		\brief consistency checks to validate that the class is set up properly for a simulation, should be called before an ensemble as in simulate()
		*/
		bool validate(double startTime, double endTime);

		/*!
		\brief selects the step size based on the propensitySum

		returns infinity if propensitySum is less than or equal to 0
		issues a warning if propensitySum is less than 0
		*/
		


		/*!
		\brief selects the index of the next reaction to fire based
		*/
		


		/*!
		\brief fire a reaction

		tested other implementations that were cleaner, but this one is faster

		\param reactionIndex the index of the reaction to fire (-1 is an error value)
		*/
		

		/*!
		\brief update the "groups" data structure

		\param affectedReactionIndex the index of the propensity that changed
		*/
		void updateGroups(int affectedReactionIndex, double oldPropensity);

	public:	
		/*!
		\brief always returns true: take a step: increment currentTime based on selectStepSize() and select and fire a reaction

		change to return type void? or return false if selectReaction returns -1 or other error occurs?

		\param reactionIndex is the number of the reaction that is firing
		*/
		bool step();

		double getCurrentTime();

		

		/*!
		\brief run an ensemble simulation with output recorded at fixed time intervals

		calls validate before ensemble
		calls initialize before each realization

		\param realizations number of simulations in the ensemble
		\param startTime the initial value of currentTime for each realization
		\param endTime the end time of each realization
		\param Output the class that handles storing the output for the simulation
		*/
		template<typename IntervalOutputType>
		void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);

	};//end SSA_ConstantTime class

    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    SSA_ConstantTime(const _populationVectorType& initialPop,
                     const _stoichiometryType& stoich,
                     const _propensitiesFunctorType& propensitiesFunctor,
                     const _dependencyGraphType& depGraph,
                     int seed) :
    initialPopulation(initialPop),
    stoichiometry(stoich),
    propensities(propensitiesFunctor),
    dependencyGraph(depGraph),
    NumberOfSpecies(initialPop.size()),
#ifdef MATRIX_STOICHIOMETRY
    NumberOfReactions(stoich.size1()),
#else
    NumberOfReactions(stoich.size()),
#endif
    currentPropensities(NumberOfReactions),
    groups(NumberOfReactions),
    initialPropensities(NumberOfReactions),
    initialGroups(NumberOfReactions)
    {
        randomGenerator.seed(seed);
        for (std::size_t i=0; i!=NumberOfReactions; ++i){
            initialPropensities[i]=propensities(i,initialPopulation);
        }
        groups.build(initialPropensities);
        initialGroups=groups;
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    seed(int seed) {
        randomGenerator.seed(seed);
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    void
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    initialize(double startTime) {
        currentTime=startTime;
        currentPopulation=initialPopulation;
        currentPropensities=initialPropensities;
        groups=initialGroups;
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    double
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    selectStepSize() {
        if (groups.getPropensitySum()<0.0) {
            //if propensitySum negative, recalculate all propensities
            groups.recalculatePropensitySum();
            //if still negative, give warning and return infinity
            if (groups.getPropensitySum()<0.0) {
                Rcpp::Rcout << "StochKit WARNING (SSA_ConstantTime::selectStepSize): propensitySum<0, returning step size=infinity\n";
                return std::numeric_limits<double>::infinity();
            }
        }
        //return randomGenerator.Exponential(1.0/groups.getPropensitySum());
        return randomGenerator.exponential(groups.getPropensitySum());
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    bool
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    validate(double startTime, double endTime) {
        if (startTime>=endTime) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): startTime not before endTime\n";
            return false;
        }
        std::size_t N=initialPopulation.size();
#ifdef MATRIX_STOICHIOMETRY
        std::size_t M=stoichiometry.size1();
#else
        std::size_t M=stoichiometry.size();
#endif
        if (N==0) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): initial population size=0\n";
            return false;
        }
        if (N!=NumberOfSpecies) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of species does not equal initial population size\n";
            return false;
        }
        if (M!=NumberOfReactions) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of reactions does not equal stoichiometry size\n";
            return false;
        }
        if (M!=propensities.size()) {
            Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): Number of reactions does not equal propensities size\n";
            return false;
        }
        //check initial populations are all non-negative
        for (std::size_t i=0; i!=NumberOfSpecies; ++i) {
            if (initialPopulation[i]<0) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): negative value detected in initial population\n";
                return false;
            }
        }
        
        //check that propensities, evaluated with initial population, are all non-negative
        for (std::size_t i=0; i!=NumberOfReactions; ++i) {
            if (propensities(i,initialPopulation)<0.0) {
                Rcpp::Rcout << "StochKit ERROR (SSA_Direct::validate): negative propensity detected based on initial population\n";
                return false;
            }
        }
        
        return true;
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    int
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    selectReaction() {
        return groups.selectReaction(randomGenerator);
    }//end selectReaction
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    bool
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    fireReaction(int reactionIndex) {
        previousReactionIndex=reactionIndex;
        
        if (reactionIndex==-1) {
            return false;
        }
        
#ifdef MATRIX_STOICHIOMETRY
        matrixrow dX(stoichiometry, reactionIndex);
        typename matrixrow::iterator it;
        for(it=dX.begin();it!=dX.end();it++) {
            currentPopulation[it.index()]+=*it;
        }
#else
        currentPopulation+=stoichiometry[reactionIndex];
#endif
        
        
        int affectedReactionIndex;
        double oldPropensity;
        
        //update affected reactions
        std::size_t numAffectedReactions=dependencyGraph[previousReactionIndex].size();
        for (std::size_t i=0; i!=numAffectedReactions; ++i) {
            affectedReactionIndex=dependencyGraph[previousReactionIndex][i];
            
            oldPropensity=currentPropensities[affectedReactionIndex];
            currentPropensities[affectedReactionIndex]=propensities(affectedReactionIndex,currentPopulation);
            
            //now we need to update groups data structures
            groups.update(affectedReactionIndex, oldPropensity, currentPropensities[affectedReactionIndex]);
        }
        
        return true;
    }//end fireReaction
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    bool
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    step() {
        currentTime+=selectStepSize();
        return fireReaction(selectReaction());
    }
    
    template<typename _populationVectorType,
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    template<typename IntervalOutputType>
    void
    SSA_ConstantTime<_populationVectorType,
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {
        
        if (doValidate) {
            if (!validate(startTime,endTime)) {
                Rcpp::Rcout << "StochKit ERROR (SSA_ConstantTime::simulate): validate() failed, simulation aborted\n";
                Rcpp::stop("Fatal error encountered, terminating StochKit2R");
            }
        }
        
        if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
            Rcpp::Rcout << "StochKit ERROR (SSA_ConstantTime::simulate): initialization of output object failed, simulation aborted\n";
            Rcpp::stop("Fatal error encountered, terminating StochKit2R");
        }
        
        std::vector<double> outputTimes = output.getOutputTimes();
        std::size_t totalIntervals=outputTimes.size();
        
        std::size_t currentInterval;
        
        for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
            initialize(startTime);
            currentInterval=0;
            
            currentTime+=selectStepSize();
            while (currentTime<endTime) {
                
                while (currentInterval<totalIntervals && currentTime>=outputTimes[currentInterval]){
                    output.record(currentRealization,currentInterval,currentPopulation);
                    currentInterval++;
                }
                
                fireReaction(selectReaction());
                currentTime+=selectStepSize();
            }
            while (currentInterval<totalIntervals && currentTime >=outputTimes[currentInterval]){
                output.record(currentRealization,currentInterval,currentPopulation);
                currentInterval++;
            }	
        }
    }//end simulate	
    
    
    template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    double
    SSA_ConstantTime<_populationVectorType, 
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    getCurrentTime(){
        return currentTime;
    }
    
    
    template<typename _populationVectorType, 
    typename _stoichiometryType,
    typename _propensitiesFunctorType,
    typename _dependencyGraphType>
    inline
    bool
    SSA_ConstantTime<_populationVectorType, 
    _stoichiometryType,
    _propensitiesFunctorType,
    _dependencyGraphType>::
    setCurrentTime(double newCurrentTime){
        currentTime = newCurrentTime;
        return true;
    }

}

#endif
