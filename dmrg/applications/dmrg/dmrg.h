/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_H
#define APP_DMRG_H

#include "dmrg/utils/DmrgParameters2.h"

#include "types/dense_matrix/dense_matrix.h"
#ifdef USE_MTM
#include "types/mt_matrix/mt_matrix.h"
typedef maquis::types::mt_matrix<double> Matrix;
#else

#ifndef USE_COMPLEX
typedef maquis::types::dense_matrix<double> Matrix;
#else
typedef maquis::types::dense_matrix<std::complex<double> > Matrix;
#endif

#endif

#include "dmrg/block_matrix/symmetry.h"


// def. of run functions
template <class SymmGroup>
void run_dmrg(DmrgParameters & parms, ModelParameters & model);

#endif
