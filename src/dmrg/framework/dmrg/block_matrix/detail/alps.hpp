/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_HPP
#define MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_HPP

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <alps/numeric/diagonal_matrix.hpp>
#include "dmrg/block_matrix/detail/alps_detail.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template <typename T, typename MemoryBlock> 
    struct transpose_view< alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> > type; }; 

    template <typename T> 
    struct transpose_view< alps::numeric::diagonal_matrix<T> > { typedef alps::numeric::diagonal_matrix<T> type; };

} }

#endif
