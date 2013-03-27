#ifndef MATRICES_H
#define MATRICES_H

#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
#endif

#include <complex>
#include <vector>

// serial matrix
#include "dmrg/block_matrix/detail/alps.hpp"

#ifdef USE_MTM
#include "types/mt_matrix/mt_matrix.h"
#include "types/mt_matrix/algorithms.hpp"
#endif

#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#include "dmrg/models/alps/model.hpp"
#endif

#include <fstream>
#include <boost/tokenizer.hpp>

// BLAS matrix
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

// MT matrix
#ifdef USE_MTM
typedef alps::numeric::mt_matrix<double> mtmatrix;
typedef alps::numeric::mt_matrix<std::complex<double> > cmtmatrix;
#endif

// parallel matrix
#ifdef USE_AMBIENT
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > pmatrix;
typedef ambient::numeric::tiles<ambient::numeric::matrix<std::complex<double> > > cpmatrix;
#endif

#endif
