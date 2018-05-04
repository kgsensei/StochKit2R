/******************************************************************************
 *  FILE:    DenseVectorSubset.h
 */

#ifndef _DENSE_VECTOR_SUBSET_H_
#define _DENSE_VECTOR_SUBSET_H_

#include <vector>
#include <fstream>

namespace STOCHKIT
{
 template<typename _denseVectorType>
 class DenseVectorSubset
 {

  public:

  DenseVectorSubset(): keepAll(true)
  {}

  DenseVectorSubset(std::vector<std::size_t> subsetIndices):keepAll(false), subsetIndices(subsetIndices)
  {}

  void setSubsetIndices(std::vector<std::size_t> subsetIndices) {
    keepAll=false;
    this->subsetIndices=subsetIndices;
  }

  void setKeepAll() {
    keepAll=true;
  }

  bool keepAll;
  std::vector<std::size_t> subsetIndices;
  

  _denseVectorType getSubset(const _denseVectorType completeVector){
    if (keepAll) {
      return completeVector;
    }
    else {
      _denseVectorType subset(subsetIndices.size());
      for (std::size_t i=0; i!=subsetIndices.size(); ++i) {
	subset[i]=completeVector[subsetIndices[i]];
      }
      return subset;
    }
  }

  std::vector<std::size_t> getSubsetIndices() {
    //will be empty if keepAll=true
    return subsetIndices;
  }
  bool getKeepAll() {
    return keepAll;
  }
  
 };
}

#endif
  
