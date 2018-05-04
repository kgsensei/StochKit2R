/*!
	\brief the explicit tau-leaping with dynamic step size selection that switches to SSA when step size is small
*/

#ifndef _TAU_LEAPING_EXPLICIT_ADAPTIVE_H_
#define _TAU_LEAPING_EXPLICIT_ADAPTIVE_H_

#include <iostream>
#include <vector>
#include <list>
#include <limits>
#include <algorithm>
#include <cmath>
#include "Random.h"
#include "StandardDriverTypes.h"
#include "SSA_DirectMatrixStoichiometry.h"
#include <boost/numeric/ublas/operation.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <Rcpp.h>

/*! 
	\param _denseVectorType the population vector type, should be dense
	\param _matrixType
	\param _propensitiesFunctorType functor takes reaction index and _denseVectorType population
           and returns the propensity for that reaction
	\param _dependencyGraphType
*/
namespace STOCHKIT
{
 template<typename _denseVectorType, 
	typename _matrixType,
	typename _propensitiesFunctorType,
	typename _dependencyGraphType>
 class TauLeapingExplicitAdaptive : public SSA_DirectMatrixStoichiometry<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType>
 {	
 public:
	//! Constructor
	TauLeapingExplicitAdaptive(const _denseVectorType& initialPop,
		   const _matrixType& stoich,
		   const _propensitiesFunctorType& propensitiesFunctor,
		   const _dependencyGraphType& depGraph,
		   int seed=time(NULL));


	//! compiler-generated copy constructor OK
	//! compiler-generated assignment operator OK

	//! destructor
	virtual ~TauLeapingExplicitAdaptive() {
	}

	void setSSASteps(std::size_t ssaSteps) {
		SSASteps=ssaSteps;
	}

	void setEpsilon(double epsilon);

	void setThreshold(std::size_t threshold) {
		this->threshold=threshold;
	}

	/*!
		\brief run an ensemble simulation with output recorded at fixed time intervals
			   
		output must have a conforming initialize, getOutputTimes, and record method
		outputTimes should be set in output prior to calling simulate
		if doValidate=true (the default) calls validate before ensemble
		calls initialize before each realization
		
		\param realizations number of simulations in the ensemble
		\param startTime the initial value of currentTime for each realization
		\param endTime the end time of each realization
		\param Output the class that handles storing the output for the simulation
	*/
	template<typename IntervalOutputType>
	void simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate=true);

 protected:
	typedef SSA_DirectMatrixStoichiometry<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType> SSA;
	using SSA::currentTime;
	using SSA::currentPopulation;
	using SSA::currentPropensities;
	using SSA::calculateAllPropensities;
	using SSA::NumberOfReactions;
	using SSA::NumberOfSpecies;
	using SSA::stoichiometry;
	using SSA::propensities;
	using SSA::randomGenerator;
	using SSA::initialPopulation;
//	using SSA::initialize;//for neg
	_denseVectorType previousReactionCounts;
	std::size_t criticalThreshold;//for neg
	double criticalPropensitySum;//for neg
	bool *tagList;//for neg
	std::vector<std::vector<std::size_t> >speciesToReaction;//for neg
	std::vector<std::vector<std::size_t> >reactionToSpecies;//for neg
	_denseVectorType affectedReactions;//for neg
	_denseVectorType affectedSpecies;//for neg

    typedef StandardDriverTypes::matrixStoichiometryRow matrixrow;
	
	std::size_t threshold;//set to zero to prevent switching to ssa
	//! the number of ssa steps taken when switching from tauleaping to ssa.
	std::size_t SSASteps;
	//! Epsilon used by tau-leaping to determine the tau
	double epsilon;

	/*! \brief the G vector described in "Efficient step size selection for the tau-leaping simulation method" */
	//_denseVectorType G;
	/* use instead of a vector--just set it to 2.0 which is conservative
		except if species can dimerize with itself and has small population */
	double g;
	//! squared elements stoichiometric matrix
	_matrixType squaredVj;
	_denseVectorType mu;
	_denseVectorType sigmaSquared;
	void prepare();
	void initialize(double startTime);
	void selectTau(double &noncriticalStepsize, double &criticalStepsize);
	
	//should this return (a reference to) the vector of reaction counts?
	int selectReactions(double leapSize, bool runCritical);

	//should this take (a reference to) the vector of reaction counts?
	bool fireReactions(int criticalIndex);
	
	//methods for critical reactions
//	double critical_selectStepSize();
//	int critical_selectReaction();
	bool critical_fireReaction(int reactionIndex);
	void critical_rollBack(int reactionIndex);
	void updateTagLists();

//	delete pure product
	std::vector<int> trimed_list;
	std::size_t NumberOfReactants;

 private:
	//! default constructor not implemented
	TauLeapingExplicitAdaptive();
	std::list<std::size_t> criticalSpecies;//for neg
	std::list<std::size_t> noncriticalSpecies;//for neg
	//! the minimum number of reactions one tauleaping step should have
	void delete_product();
 };//end TauLeapingExplicitAdaptive class

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::
TauLeapingExplicitAdaptive(const _denseVectorType& initialPop,
                           const _matrixType& stoich,
                           const _propensitiesFunctorType& propensitiesFunctor,
                           const _dependencyGraphType& depGraph,
                           int seed) : SSA_DirectMatrixStoichiometry<_denseVectorType, _matrixType, _propensitiesFunctorType, _dependencyGraphType>(initialPop, stoich, propensitiesFunctor, depGraph, seed),
previousReactionCounts(NumberOfReactions),
criticalThreshold(5),
affectedReactions(NumberOfSpecies),
affectedSpecies(NumberOfReactions),
threshold(10),
SSASteps(100),
epsilon(0.03),
g(3.0),//could improve
squaredVj(stoich),
mu(NumberOfSpecies),
sigmaSquared(NumberOfSpecies)
{
    squaredVj=element_prod(stoichiometry, stoichiometry);
    //should have a debug check here that ensures _denseVectorType::value_type and _matrixType::value_type are both double
    //for neg
    tagList=new bool[NumberOfReactions];
    std::vector<std::size_t> temp;
    std::size_t i, j;
    for(i=0; i<NumberOfReactions; i++)
        reactionToSpecies.push_back(temp);
    for(i=0; i<NumberOfSpecies; i++)
        speciesToReaction.push_back(temp);
    for(i=0; i<NumberOfReactions; i++)
    {
        for(j=0; j<NumberOfSpecies; j++)
        {
            if(stoichiometry(i,j)<0)
            {
                reactionToSpecies[i].push_back(j);
                speciesToReaction[j].push_back(i);
            }
        }
    }
    for(i=0; i<NumberOfReactions; i++)
        affectedSpecies[i]=reactionToSpecies[i].size();
    for(i=0; i<NumberOfSpecies; i++)
        affectedReactions[i]=speciesToReaction[i].size();
}


template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::delete_product()//for neg
{
    trimed_list.clear();
    for(std::size_t i=0; i<NumberOfSpecies; i++)
    {
        if(affectedReactions[i]!=0)
            trimed_list.push_back(i);
    }
    NumberOfReactants=trimed_list.size();
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::initialize(double startTime)//for neg
{
    std::size_t i;
    SSA::initialize(startTime);
    criticalSpecies.clear();
    noncriticalSpecies.clear();
    for(i=0; i<NumberOfReactants; i++)
        noncriticalSpecies.push_back(trimed_list[i]);
    for(i=0; i<NumberOfReactions; i++)
        tagList[i]=false;
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::prepare()//for neg
{
    delete_product();
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::setEpsilon(double epsilon)
{
    if (epsilon<=0.0 || epsilon>=1.0) {
        Rcpp::Rcout << "StochKit ERROR (TauLeapingExplicitAdaptive::setEpsilon): epsilon must be greater than 0.0 and less than 1.0.  Terminating simulation.\n";
        Rcpp::stop("Fatal error encountered. Terminating");
    }
    this->epsilon=epsilon;
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::selectTau(double &noncriticalStepsize, double &criticalStepsize)
{
    //for critical reactions
    criticalPropensitySum=0;
    for (std::size_t i=0; i!=NumberOfReactions; ++i) {
        if(tagList[i]==true)
            criticalPropensitySum+=currentPropensities[i];
    }
    if(criticalPropensitySum!=0)
        criticalStepsize=randomGenerator.exponential(criticalPropensitySum);
    else
        criticalStepsize=-1;
    
    //for noncritical reactions
    //updateMuAndSigmaSquared();
    axpy_prod(currentPropensities, stoichiometry, mu, true);
    axpy_prod(currentPropensities, squaredVj, sigmaSquared, true);
    
    noncriticalStepsize=std::numeric_limits<double>::max();
    double numerator, temp1, temp2;
    std::vector<int>::iterator Iterator;
    for(Iterator=trimed_list.begin(); Iterator!=trimed_list.end(); Iterator++)
    {
        numerator=epsilon*currentPopulation[*Iterator]/g;//temp--should use G
        numerator=std::max(numerator,1.0);
        temp1=numerator/fabs(mu(*Iterator));
        temp2=numerator*numerator/sigmaSquared(*Iterator);
        noncriticalStepsize=std::min(noncriticalStepsize,temp1);
        noncriticalStepsize=std::min(noncriticalStepsize,temp2);
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
int
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::selectReactions(double leapSize, bool runCritical)
{
    std::size_t i;
    if(runCritical)//handle both critical and noncritical reactions
    {
        //for critical
        std::size_t previousReactionIndex=-1;
        //generate a uniform random number between (0,propensitySum)
        double r=0;
        while (r==0) {
            r=randomGenerator.continuousUniform(0,1)*criticalPropensitySum;
        }
        double jsum=0;
        
        for (i=0; i!=NumberOfReactions; ++i) {
            //for critical
            if(tagList[i]==true)
            {
                previousReactionCounts[i]=0;//not noncritical so set 0
                if(jsum<r)
                {
                    previousReactionIndex=i;
                    jsum+=currentPropensities[i];
                }
            }
            //for noncritical
            else
                previousReactionCounts[i]=randomGenerator.poisson(leapSize*currentPropensities[i]);
        }
        return previousReactionIndex;
    }
    else//only handle noncritical reactions
    {
        for (i=0; i!=NumberOfReactions; ++i) {
            if(tagList[i]==true)
                previousReactionCounts[i]=0;
            else
                previousReactionCounts[i]=randomGenerator.poisson(leapSize*currentPropensities[i]);
        }
        return -1;
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::fireReactions(int criticalIndex)
{
    _denseVectorType populationChange(NumberOfSpecies);
    axpy_prod(previousReactionCounts,stoichiometry,populationChange,true);
    currentPopulation+=populationChange;
    if(criticalIndex!=-1)
        critical_fireReaction(criticalIndex);
    
    //check for negative populations
    bool negativePopulation=false;
    for (std::size_t i=0; i!=NumberOfSpecies; ++i) {
        if (currentPopulation[i]<0.0) {
            negativePopulation=true;
            break;
        }
    }
    if (negativePopulation) {
        currentPopulation-=populationChange;
        if(criticalIndex!=-1)
            critical_rollBack(criticalIndex);
        return false;
    }
    else {
        calculateAllPropensities();
        return true;
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
bool
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::
critical_fireReaction(int reactionIndex) {
    if (reactionIndex==-1) {
        return false;
    }
    else {
        //if not -1, assumes valid reactionIndex
        matrixrow dX(stoichiometry, reactionIndex);
        typename matrixrow::iterator it;
        for(it=dX.begin();it!=dX.end();it++) {
            currentPopulation[it.index()]+=*it;
        }
        return true;
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::
critical_rollBack(int reactionIndex) {
    //if not -1, assumes valid reactionIndex
    matrixrow dX(stoichiometry, reactionIndex);
    typename matrixrow::iterator it;
    for(it=dX.begin();it!=dX.end();it++) {
        currentPopulation[it.index()]-=*it;
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::
updateTagLists()
{
    std::size_t j, k;
    bool changeTag;
    std::list<std::size_t>::iterator Iterator;
    
    //from critical species to noncritical species
    for(Iterator=criticalSpecies.begin(); Iterator!=criticalSpecies.end();)
    {
        if(currentPopulation[*Iterator]>criticalThreshold)
        {
            noncriticalSpecies.push_back(*Iterator);//add to noncritical species list
            for(j=0; j<affectedReactions[*Iterator]; j++)//modify tag for each reaction that take this species as a reactant if neccessary
            {
                changeTag=true;
                for(k=0; k<affectedSpecies[speciesToReaction[*Iterator][j]]; k++)//check if other reactants are critical for this reaction
                {
                    if(currentPopulation[reactionToSpecies[speciesToReaction[*Iterator][j]][k]]<=criticalThreshold)
                    {
                        changeTag=false;
                        break;
                    }
                }
                if(changeTag)
                    tagList[speciesToReaction[*Iterator][j]]=false;
            }
            Iterator=criticalSpecies.erase(Iterator);//erase from critical list
        }
        else
            Iterator++;
    }
    
    //from noncritical species to critical species
    for(Iterator=noncriticalSpecies.begin(); Iterator!=noncriticalSpecies.end();)
    {
        if(currentPopulation[*Iterator]<=criticalThreshold)
        {
            criticalSpecies.push_back(*Iterator);//add to critical species list
            for(j=0; j<affectedReactions[*Iterator]; j++)//modify tag for each reaction that take this species as a reactant
                tagList[speciesToReaction[*Iterator][j]]=true;
            Iterator=noncriticalSpecies.erase(Iterator);
        }
        else
            Iterator++;
    }
}

template<typename _denseVectorType,
typename _matrixType,
typename _propensitiesFunctorType,
typename _dependencyGraphType>
template<typename IntervalOutputType>
void
TauLeapingExplicitAdaptive<_denseVectorType,
_matrixType,
_propensitiesFunctorType,
_dependencyGraphType>::simulate(std::size_t realizations, double startTime, double endTime, IntervalOutputType& output, bool doValidate) {
    
    //validate();
    if (doValidate) {
        if (!SSA::validate(startTime,endTime)) {
            Rcpp::Rcout << "StochKit ERROR (TauLeapingExplicitAdaptive::simulate): validate() failed, simulation aborted\n";
            Rcpp::stop("Fatal error encountered. Terminating");
        }
    }
    
    if (!output.initialize(realizations,startTime,endTime,initialPopulation)) {
        Rcpp::Rcout << "StochKit ERROR (): initialization of output object failed, simulation aborted\n";
        Rcpp::stop("Fatal error encountered. Terminating");
    }
    
    //delete the product for tau selection
    prepare();
    
    std::vector<double> outputTimes=output.getOutputTimes();
    std::size_t totalIntervals=outputTimes.size();
    std::size_t currentInterval;
    std::size_t failedLeaps;
    double tau;
    double nextTime;
    double reactionsLastLeap;
    std::size_t ssaStepsTaken;
    double criticalStepsize;//for neg
    double noncriticalStepsize;//for neg
    bool runCritical;//for neg
    int criticalIndex;//for negi
    
    for (std::size_t currentRealization=0; currentRealization!=realizations; ++currentRealization) {
        initialize(startTime);
        currentInterval=0;
        failedLeaps=0;
        ssaStepsTaken=0;
        reactionsLastLeap=std::numeric_limits<double>::max();
        
        while (currentTime<endTime) {
            while ((currentTime>=outputTimes[currentInterval]) && currentInterval<totalIntervals) {
                output.record(currentRealization, currentInterval, currentPopulation);
                ++currentInterval;
            }
            
            if (reactionsLastLeap<threshold) {
                //do SSA
                failedLeaps=0;//not taking leaps
                if (ssaStepsTaken<SSASteps) {
                    if (ssaStepsTaken==0) {
                        currentTime+=SSA::selectStepSize();
                        ++ssaStepsTaken;
                    }
                    else {
                        SSA::fireReaction(SSA::selectReaction());
                        currentTime+=SSA::selectStepSize();
                        ++ssaStepsTaken;
                    }
                }
                else { // we've taken enough SSA steps, try to go back to tau-leaping
                    SSA::fireReaction(SSA::selectReaction());
                    reactionsLastLeap=std::numeric_limits<double>::max();
                    ssaStepsTaken=0;
                }
            }
            else {
                //do tau-leaping
                updateTagLists();//for neg
                selectTau(noncriticalStepsize, criticalStepsize);//for neg
                if(noncriticalStepsize<criticalStepsize||criticalStepsize==-1)//for neg
                {
                    tau=noncriticalStepsize;
                    runCritical=false;
                }
                else
                {
                    tau=criticalStepsize;
                    runCritical=true;
                }
                if (currentTime+tau>outputTimes[currentInterval]) {
                    //don't leap past next output time
                    tau=outputTimes[currentInterval]-currentTime;
                    nextTime=outputTimes[currentInterval];
                    runCritical=false;//for neg
                }
                else {
                    nextTime=currentTime+tau;
                }
                criticalIndex=selectReactions(tau, runCritical);
                reactionsLastLeap=norm_1(previousReactionCounts);
                if(runCritical)//for neg
                    reactionsLastLeap++;
                if (fireReactions(criticalIndex)) {
                    currentTime=nextTime;
                    failedLeaps=0;
                }
                else {
                    ++failedLeaps;
                    if (failedLeaps==3) {
                        Rcpp::Rcout << "StochKit WARNING (TauLeapingExplicitAdaptive::simulate): rejected three or more consecutive leaps, consider reducing epsilon.\n";
                    }
                    if (failedLeaps==10) {
                        Rcpp::Rcout << "StochKit ERROR (TauLeapingExplicitAdaptive::simulate): rejected ten consecutive leaps, terminating simulation.\n";
                        Rcpp::stop("Fatal error encountered. Terminating");
                    }
                }
            }
        }
        
        //modified for visual studio
        while(currentTime>=outputTimes[currentInterval])
        {
            output.record(currentRealization, currentInterval, currentPopulation);
            ++currentInterval;
            if(currentInterval>=totalIntervals)
                break;
        }
    }
}//end simulate

}//end namespace

#endif
