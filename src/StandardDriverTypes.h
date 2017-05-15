/******************************************************************************
 */

#ifndef _STANDARD_TYPES_H_
#define _STANDARD_TYPES_H_
#include "boost/numeric/ublas/vector_sparse.hpp"
#include "boost/numeric/ublas/vector.hpp"
#include "boost/numeric/ublas/matrix_sparse.hpp"
#include "boost/numeric/ublas/matrix_proxy.hpp"
#include <vector>
#include "CustomPropensitySet.h"
//#include "CustomPropensitySet_Events.h"
#include "output/StandardDriverOutput.h"

namespace STOCHKIT
{
 class StandardDriverTypes
 {
	
 public:

  typedef boost::numeric::ublas::vector<double> populationType;

  typedef boost::numeric::ublas::mapped_vector<double> sparse_double_ublas_vec;

  typedef boost::numeric::ublas::vector<double> dense_double_ublas_vec;

    typedef boost::numeric::ublas::compressed_matrix<double> matrixStoichiometryType;
    typedef boost::numeric::ublas::matrix_row<matrixStoichiometryType> matrixStoichiometryRow;
    typedef std::vector<sparse_double_ublas_vec> sparseStoichiometryType;
    typedef std::vector<dense_double_ublas_vec> denseStoichiometryType;//for "small" drivers

  typedef std::vector<std::vector<std::size_t> > graphType;
  typedef CustomPropensitySet<populationType> propensitiesType;
//  typedef CustomPropensitySet_Events<populationType> eventEnabledPropensitiesType;
  typedef StandardDriverOutput<populationType> outputType;

 };
}

#endif
