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
#include "dmrg/block_matrix/detail/one_matrix.hpp"
#include "utils/traits.hpp"

namespace maquis { namespace traits {

    template <typename T, typename MemoryBlock> 
    struct transpose_view< alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::transpose_view<alps::numeric::matrix<T, MemoryBlock> > type; }; 

    template <typename T> 
    struct transpose_view< alps::numeric::diagonal_matrix<T> > { typedef alps::numeric::diagonal_matrix<T> type; };

} }

namespace alps { namespace numeric {

    template<class Matrix> struct associated_one_matrix { };
    template<class Matrix> struct associated_dense_matrix { };

    template<typename T, typename MemoryBlock>
    struct associated_one_matrix<alps::numeric::matrix<T, MemoryBlock> > { typedef maquis::dmrg::one_matrix<T> type; };

    template<typename T, class MemoryBlock>
    struct associated_dense_matrix<alps::numeric::matrix<T, MemoryBlock> > { typedef alps::numeric::matrix<T, MemoryBlock> type; };

    template<typename T>
    struct associated_dense_matrix<maquis::dmrg::one_matrix<T> > { typedef alps::numeric::matrix<T> type; };

} }

#endif
