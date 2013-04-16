#ifndef MATRICES_H
#define MATRICES_H

#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
#endif

#include <complex>
#include <vector>

// serial matrix
#include "dmrg/block_matrix/detail/alps.hpp"

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#include "dmrg/models/alps/model.hpp"
#endif

#include <fstream>
#include <boost/tokenizer.hpp>

// BLAS matrix
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

// parallel matrix
#ifdef USE_AMBIENT
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > pmatrix;
typedef ambient::numeric::tiles<ambient::numeric::matrix<std::complex<double> > > cpmatrix;
#endif

#endif
