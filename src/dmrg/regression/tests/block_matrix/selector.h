/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TEST_DMRG_H
#define TEST_DMRG_H

#ifdef USE_COMPLEX
    #define dmrg_value_type std::complex<double>
#else
    #define dmrg_value_type double
#endif

#ifdef PDENSE
    #include "types/p_dense_matrix/p_dense_matrix.h"
    #include "types/p_dense_matrix/algorithms.hpp"
    #include "dmrg/block_matrix/detail/ambient_matrix_kernels.hpp"
    typedef maquis::types::p_dense_matrix<dmrg_value_type> Matrix;
#endif

#ifdef DENSE
    #include "dmrg/block_matrix/detail/alps_matrix.hpp"
    typedef alps::numeric::matrix<dmrg_value_type> Matrix;
#endif

#undef dmrg_value_type

#endif
