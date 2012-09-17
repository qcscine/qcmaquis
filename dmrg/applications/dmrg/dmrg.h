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
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> > matrix;
#elif defined USE_MTM
    #include "types/mt_matrix/mt_matrix.h"
    typedef alps::numeric::mt_matrix<dmrg_value_type> matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<dmrg_value_type> matrix;
#endif
#undef dmrg_value_type

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/block_matrix/symmetry.h"

// def. of run functions
template <class SymmGroup>
void run_dmrg(DmrgParameters & parms, ModelParameters & model);

#endif
