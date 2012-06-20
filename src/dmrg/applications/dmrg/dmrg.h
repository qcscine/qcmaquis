/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_DMRG_H
#define APP_DMRG_H

#ifdef USE_COMPLEX
    #define dmrg_value_type std::complex<double>
#else
    #define dmrg_value_type double
#endif
#if defined USE_AMBIENT
    #include "types/p_dense_matrix/p_dense_matrix.h"
    typedef maquis::types::p_dense_matrix<dmrg_value_type> matrix;
#elif defined USE_MTM
    #include "types/mt_matrix/mt_matrix.h"
    typedef alps::numeric::mt_matrix<dmrg_value_type> matrix;
#else
    #include "alps/numeric/matrix/matrix.hpp"
    #include "alps/numeric/matrix/matrix_interface.hpp"
    #include "alps/numeric/matrix/resizable_matrix_interface.hpp"
    #include "alps/numeric/matrix/algorithms.hpp"
    #include "alps/numeric/matrix/matrix_blas.hpp"
    typedef alps::numeric::matrix<dmrg_value_type> matrix;
#endif
#undef dmrg_value_type

#include "dmrg/utils/DmrgParameters2.h"
#include "alps/numeric/matrix/matrix.hpp"
#include "dmrg/block_matrix/symmetry.h"

// def. of run functions
template <class SymmGroup>
void run_dmrg(DmrgParameters & parms, ModelParameters & model);

#endif
