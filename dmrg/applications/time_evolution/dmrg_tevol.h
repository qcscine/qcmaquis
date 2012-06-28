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
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::matrix<double> matrix;
    typedef ambient::numeric::matrix<std::complex<double> > cmatrix;
#elif defined USE_MTM
#include "types/mt_matrix/mt_matrix.h"
    typedef alps::numeric::mt_matrix<double> matrix;
    typedef alps::numeric::mt_matrix<std::complex<double> > cmatrix;
#else
    #include "dmrg/block_matrix/detail/alps_matrix.hpp"
    typedef alps::numeric::matrix<double> matrix;
    typedef alps::numeric::matrix<std::complex<double> > cmatrix;
#endif

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/block_matrix/detail/alps_matrix.hpp"
#include "dmrg/block_matrix/symmetry.h"

// def. of run functions
template <class Matrix, class SymmGroup>
void run_tevol(DmrgParameters & parms, ModelParameters & model);

#endif
