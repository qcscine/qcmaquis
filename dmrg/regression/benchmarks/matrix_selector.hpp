/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BENCHMARKS_MATRIX_SELECTOR_HPP
#define BENCHMARKS_MATRIX_SELECTOR_HPP

#include <vector>

// serial matrix
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> amatrix;

// parallel matrix
#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > pmatrix;
#endif

#ifdef USE_AMBIENT
typedef pmatrix matrix;
#else
typedef amatrix matrix;
#endif

#endif
