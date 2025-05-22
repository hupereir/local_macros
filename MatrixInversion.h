#ifndef MatrixInversion_h
#define MatrixInversion_h

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>

namespace ublas = boost::numeric::ublas;

/* Matrix inversion routine. Uses lu_factorize and lu_substitute in uBLAS to invert a matrix */
template<class T>
bool invert (const boost::numeric::ublas::matrix<T>& input, boost::numeric::ublas::matrix<T>& inverse)
{
  using namespace boost::numeric::ublas;
  typedef permutation_matrix<std::size_t> pmatrix;

  // create a working copy of the input
  matrix<T> A(input);

  // create a permutation matrix for the LU-factorization
  pmatrix pm(A.size1());

  // perform LU-factorization
  auto res = lu_factorize(A,pm);
  if( res != 0 ) return false;

  // create identity matrix of "inverse"
  inverse = ublas::identity_matrix<T>(input.size1());

  // backsubstitute to get the inverse
  lu_substitute(A, pm, inverse);

  return true;
}

#endif
