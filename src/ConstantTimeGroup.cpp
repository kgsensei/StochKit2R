#include "ConstantTimeGroup.h"

namespace STOCHKIT
{

ConstantTimeGroup::ConstantTimeGroup(int groupExponent): groupSum(0), groupExponent(groupExponent), maxPropensity(pow(2.0,(double)groupExponent)), minPropensity(pow(2.0,(double)(groupExponent-1)))
{
}

double ConstantTimeGroup::getGroupSum() {
  return groupSum;
}

int ConstantTimeGroup::getGroupExponent() {
  return groupExponent;
}

int ConstantTimeGroup::calculateGroupExponent(double propensityValue) {
  int exponent=(int)ceil(log2(propensityValue));
  return exponent;
}

int ConstantTimeGroup::selectReactionIndex(Random &randomGenerator) {
  int tentativePropensityValuesIndex=-1;

  double r;
  
  if (reactionIndexes.size()==1) { //if there is only one reaction in the group, choose it
    return reactionIndexes[0];
  }
  else {
    while (true) {
      //generate an integer index between 0 and size of group-1 as tentative reaction index 
      tentativePropensityValuesIndex=(int)(randomGenerator.continuousUniform(0,1)*reactionIndexes.size());

      //generate a continuous number within the group's propensity range
      r=randomGenerator.continuousUniform(0,1)*maxPropensity;
      //if r is less than the propensity of tentativeReaction propensity then accept the reaction
      if (r<=propensityValues[tentativePropensityValuesIndex]) {
	return (int)reactionIndexes[tentativePropensityValuesIndex];
      }
    }
  }
  
  //should never get here
  return -1;
}

void ConstantTimeGroup::insert(double propensityValue, std::size_t reactionIndex, std::vector<int> &groupsWithinGroupIndexes) {

  groupSum+=propensityValue;
  propensityValues.push_back(propensityValue);
  reactionIndexes.push_back(reactionIndex);

  groupsWithinGroupIndexes[reactionIndex]=reactionIndexes.size()-1;
}

void ConstantTimeGroup::remove(std::size_t reactionIndex, std::size_t withinGroupIndex, std::vector<int> &groupsWithinGroupIndexes) {

  groupSum-=propensityValues[withinGroupIndex];
  
  groupsWithinGroupIndexes[reactionIndexes.back()]=withinGroupIndex;
  groupsWithinGroupIndexes[reactionIndex]=-1;

  propensityValues[withinGroupIndex]=propensityValues.back();
  reactionIndexes[withinGroupIndex]=reactionIndexes.back();
  reactionIndexes.pop_back();
  propensityValues.pop_back();

  if (propensityValues.size()==0) {
    groupSum=0;
  }
}

void ConstantTimeGroup::update(double newPropensityValue, std::size_t reactionIndex, std::size_t withinGroupIndex, std::vector<int> &groupsWithinGroupIndexes) {

  groupSum+=newPropensityValue-propensityValues[withinGroupIndex];
  propensityValues[withinGroupIndex]=newPropensityValue;
}


}
