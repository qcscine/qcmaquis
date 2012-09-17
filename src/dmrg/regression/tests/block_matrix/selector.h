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
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<dmrg_value_type> > Matrix;
#endif

#ifdef DENSE
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<dmrg_value_type> Matrix;
#endif

#undef dmrg_value_type

#endif
