#include "ConstantTimeGroupCollection.h"
#include <Rcpp.h>

namespace STOCHKIT
{

ConstantTimeGroupCollection::ConstantTimeGroupCollection(std::size_t numberOfReactions): propensitySum(0), withinGroupIndexes(numberOfReactions,-1)
{
  maxGroupExponent=std::numeric_limits<int>::min();
  minGroupExponent=std::numeric_limits<int>::max();
  if (numberOfReactions==0) {
      Rcpp::Rcout << "StochKit ERROR (ConstantTimeGroupCollection): must have at least one reaction. Terminating.\n";
      Rcpp::stop("Fatal error occurred, terminating StochKit2R");
  }
}

double ConstantTimeGroupCollection::getPropensitySum() {
  return propensitySum;
}

void ConstantTimeGroupCollection::build(std::vector<double> &propensities) {
  groups.clear();

  bool oneNonZeroPropensity=false;

  for (std::size_t i=0; i!=propensities.size(); ++i) {
    if (!oneNonZeroPropensity) {
      if (propensities[i]>0.0) {
	oneNonZeroPropensity=true;
	int exponent=ConstantTimeGroup::calculateGroupExponent(propensities[i]);
	minGroupExponent=exponent;
	maxGroupExponent=exponent;
	ConstantTimeGroup newGroup(exponent);
	groups.push_front(newGroup);
	groups[0].insert(propensities[i],i,withinGroupIndexes);
	propensitySum=propensities[i];
      }
    }
    else {
      if (propensities[i]>0.0) {
	update(i,0.0,propensities[i]);
      }
    }
  }

  if (!oneNonZeroPropensity) {
      Rcpp::Rcout << "StochKit ERROR (ConstantTimeGroupCollection::build): All propensities are zero\n";
      Rcpp::stop("Fatal error occurred, terminating StochKit2R");
  }
  
}

void ConstantTimeGroupCollection::update(std::size_t reactionIndex, double oldPropensity, double newPropensity) {
  propensitySum+=newPropensity-oldPropensity;

  int newGroup=getGroup(newPropensity);
  int oldGroup=getGroup(oldPropensity);

  if (newGroup==oldGroup) {
    if (newGroup==-1) {
      //either propensity=0 (do nothing) or need to add new group
      if (newPropensity!=0.0) {
	addGroup(ConstantTimeGroup::calculateGroupExponent(newPropensity));
	newGroup=getGroup(newPropensity);//group index changed
	groups[newGroup].insert(newPropensity,reactionIndex,withinGroupIndexes);
      }
    }
    else {
      //did not change group, simple update
      groups[newGroup].update(newPropensity, reactionIndex, withinGroupIndexes[reactionIndex], withinGroupIndexes);
    }
  }
  else {//changed group
    //remove from old group
    if (oldGroup!=-1) {
      groups[oldGroup].remove(reactionIndex,withinGroupIndexes[reactionIndex],withinGroupIndexes);
    }
    //add to new group
    if (newGroup==-1) {
      if (newPropensity>0) {
	//need to add a group
	addGroup(ConstantTimeGroup::calculateGroupExponent(newPropensity));
	newGroup=getGroup(newPropensity);//group index changed
	groups[newGroup].insert(newPropensity,reactionIndex,withinGroupIndexes);
      }
    }
    else {
      groups[newGroup].insert(newPropensity, reactionIndex, withinGroupIndexes);
    }
  }

}

int ConstantTimeGroupCollection::selectReaction(Random &randomGenerator) {
  //first select a group via linear search
  int groupIndex=selectGroupIndex(randomGenerator);
  
  if (groupIndex==-1) {
    return -1;
  }
  else {
    int reactionIndex=groups[groupIndex].selectReactionIndex(randomGenerator);
    return reactionIndex;
  }
}

void ConstantTimeGroupCollection::recalculatePropensitySum() {
  propensitySum=0.0;
  for (std::size_t i=0; i!=groups.size(); ++i) {
    propensitySum+=groups[i].getGroupSum();
  }
}

int ConstantTimeGroupCollection::getGroup(double propensity) {
  if (propensity==0) {
    return -1;
  }
  else {
    int exponent=ConstantTimeGroup::calculateGroupExponent(propensity);
    if (exponent>=minGroupExponent && exponent<=maxGroupExponent) {
      return maxGroupExponent-exponent;
    }
    else {
      return -1;
    }
  }
}

void ConstantTimeGroupCollection::addGroup(int newGroupExponent) {
  while (newGroupExponent<minGroupExponent) {
    ConstantTimeGroup newGroup(--minGroupExponent);
    groups.push_back(newGroup);
  }
  while (newGroupExponent>maxGroupExponent) {
    ConstantTimeGroup newGroup(++maxGroupExponent);
    groups.push_front(newGroup);    
  }
}

int ConstantTimeGroupCollection::selectGroupIndex(Random &randomGenerator) {
  int groupIndex=-1;
  double r=0;
  while (r==0) {
	r=randomGenerator.continuousUniform(0,1)*propensitySum;
  }
  double jsum=0;
  while (jsum < r) {
    ++groupIndex;
    if (groupIndex==(int)groups.size()) {
      recalculatePropensitySum();
      return selectGroupIndex(randomGenerator);
    }
    else {
      jsum+=groups[groupIndex].getGroupSum();
    }
  }
  return groupIndex;
}

}

