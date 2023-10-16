/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SIM_MATRIX_TYPES_H
#define SIM_MATRIX_TYPES_H

#if defined USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
#include "dmrg/block_matrix/detail/alps.hpp"
#include <complex>
typedef ambient::tiles<ambient::matrix<double> > matrix;
typedef ambient::tiles<ambient::matrix< std::complex<double> > > cmatrix;
template <class V>
    using tmatrix = ambient::tiles<ambient::matrix<V> >;
#else
#include "dmrg/block_matrix/detail/alps.hpp"
#include <complex>
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;
template <class V>
    using tmatrix = alps::numeric::matrix<V>;
#endif

#endif
