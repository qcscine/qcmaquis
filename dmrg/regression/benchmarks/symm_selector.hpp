/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BENCHMARKS_SYMM_SELECTOR_HPP
#define BENCHMARKS_SYMM_SELECTOR_HPP

#include "dmrg/block_matrix/symmetry.h"

#ifdef USE_TWOU1
typedef TwoU1 grp;
#else
#ifdef USE_NONE
typedef TrivialGroup grp;
#else
typedef U1 grp;
#endif
#endif

#endif
