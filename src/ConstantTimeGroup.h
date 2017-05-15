#ifndef _CONSTANT_TIME_GROUP_H_
#define _CONSTANT_TIME_GROUP_H_

#include <iostream>
#include <sstream>
#include <vector>
#include <string>
#include <utility>
#include <limits>
#include <cmath>
#include "Random.h"

namespace STOCHKIT
{
 class ConstantTimeGroup
 {
 public:

  ConstantTimeGroup(int groupExponent);

  void insert(double propensityValue, std::size_t reactionIndex, std::vector<int> &groupsWithinGroupIndexes);
  void remove(std::size_t reactionIndex, std::size_t withinGroupIndex, std::vector<int> &groupsWithinGroupIndexes);

  void update(double newPropensityValue, std::size_t reactionIndex, std::size_t withinGroupIndex, std::vector<int> &groupsWithinGroupIndexes);

  int selectReactionIndex(Random &randomGenerator);

  double getGroupSum();
  int getGroupExponent();

  static int calculateGroupExponent(double propensityValue);


 private:
  ConstantTimeGroup();

  double groupSum;
  std::vector<double> propensityValues;
  std::vector<std::size_t> reactionIndexes;
  int groupExponent;
  double maxPropensity;
  double minPropensity;//lower limit is not included in this group

 };
}

#endif

