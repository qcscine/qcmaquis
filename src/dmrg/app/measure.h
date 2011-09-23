/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MEASURE_H
#define APP_MEASURE_H

#include "utils/DmrgParameters2.h"

#include "dense_matrix/dense_matrix.h"
#ifdef USE_MTM
#include "dense_matrix/mt_matrix.h"
typedef mt_matrix<double> Matrix;
#else
typedef blas::dense_matrix<double> Matrix;
#endif

#include "block_matrix/symmetry.h"


// def. of run functions
template <class SymmGroup>
void run_dmrg(DmrgParameters & parms, ModelParameters & model);


#endif
