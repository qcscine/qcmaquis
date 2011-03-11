/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SIZEOF_H
#define SIZEOF_H

#include "dense_matrix/dense_matrix.h"
#include "block_matrix/block_matrix.h"
#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

namespace utils
{
    using std::size_t;
    
    template<class T, class MB>
    size_t size_of(blas::dense_matrix<T, MB> const & m)
    {
        return num_rows(m)*num_columns(m)*sizeof(T);
    }
    
    template<class M, class SG>
    size_t size_of(block_matrix<M, SG> const & m)
    {
        size_t r = 0;
        for (size_t i = 0; i < m.n_blocks(); ++i)
            r += size_of(m[i]);
        return r;
    }
    
    template<class M, class SG>
    size_t size_of(MPSTensor<M, SG> const & m)
    {
        return size_of(m.data());
    }
    
    template<class M, class SG>
    size_t size_of(Boundary<M, SG> const & m)
    {
        size_t r = 0;
        for (size_t i = 0; i < m.aux_dim(); ++i)
            r += size_of(m[i]);
        return r;
    }
    
    template<class Iterator>
    size_t size_of(Iterator i1, Iterator i2)
    {
        size_t r = 0;
        for ( ; i1 != i2; ++i1)
            r += size_of(*i1);
        return r;
    }
}

#endif
