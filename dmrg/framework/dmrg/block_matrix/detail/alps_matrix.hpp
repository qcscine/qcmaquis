/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_H
#define MAQUIS_BLOCK_MATRIX_DEATAIL_ALPS_H

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include <alps/numeric/diagonal_matrix.hpp>
#include "dmrg/block_matrix/detail/alps_matrix_kernels.hpp"

#include "types/utils/traits.hpp"


namespace maquis { namespace traits {

    template <typename T, typename MemoryBlock>
    struct transpose< alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> > type; };
    template <typename T, typename MemoryBlock>
    struct transpose< alps::numeric::matrix<T, MemoryBlock> const > { typedef alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> const > type; };
    
    // MD: this is wrong, we need a workaround!
    template <typename T, typename MemoryBlock>
    struct transpose< alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> > > { typedef alps::numeric::matrix<T, MemoryBlock> type; };
    template <typename T, typename MemoryBlock>
    struct transpose< alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> const > > { typedef alps::numeric::matrix<T, MemoryBlock> const type; };

    template <typename T>
    struct transpose< alps::numeric::diagonal_matrix<T> > { typedef alps::numeric::diagonal_matrix<T> type; };

} }

#endif
