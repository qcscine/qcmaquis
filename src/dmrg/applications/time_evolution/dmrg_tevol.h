/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_TEVOL_H
#define APP_DMRG_TEVOL_H

#if defined USE_AMBIENT
    #include "types/p_dense_matrix/p_dense_matrix.h"
    typedef maquis::types::p_dense_matrix<double> matrix;
    typedef maquis::types::p_dense_matrix<std::complex<double> > cmatrix;
#elif defined USE_MTM
#include "types/mt_matrix/mt_matrix.h"
    typedef maquis::types::mt_matrix<double> matrix;
    typedef maquis::types::mt_matrix<std::complex<double> > cmatrix;
#else
    #include "types/dense_matrix/dense_matrix.h"
    #include "types/dense_matrix/matrix_interface.hpp"
    #include "types/dense_matrix/resizable_matrix_interface.hpp"
    #include "types/dense_matrix/dense_matrix_blas.hpp"
    #include "types/dense_matrix/algorithms.hpp"
    typedef maquis::types::dense_matrix<double> matrix;
    typedef maquis::types::dense_matrix<std::complex<double> > cmatrix;
#endif

#include "dmrg/utils/DmrgParameters2.h"
#include "types/dense_matrix/dense_matrix.h"
#include "dmrg/block_matrix/symmetry.h"

// def. of run functions
template <class Matrix, class SymmGroup>
void run_tevol(DmrgParameters & parms, ModelParameters & model);

#endif
