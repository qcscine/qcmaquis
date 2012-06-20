/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_MEASURE_H
#define APP_MEASURE_H

#include "dmrg/utils/DmrgParameters2.h"

#include "alps/numeric/matrix/matrix.hpp"
#include "alps/numeric/matrix/matrix_interface.hpp"
#include "alps/numeric/matrix/resizable_matrix_interface.hpp"
#include "alps/numeric/matrix/matrix_blas.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#ifdef USE_MTM
#include "types/mt_matrix/mt_matrix.h"
typedef alps::numeric::mt_matrix<double> Matrix;
#else
typedef alps::numeric::matrix<double> Matrix;
#endif

#include "dmrg/block_matrix/symmetry.h"


// def. of run functions
template <class SymmGroup>
void run_mg_meas(DmrgParameters & parms, ModelParameters & model);


#endif
