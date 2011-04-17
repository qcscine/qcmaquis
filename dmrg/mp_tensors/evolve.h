/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef EVOLVE_H
#define EVOLVE_H

#include "utils/zout.hpp"

#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

#include "mp_tensors/reshapes.h"
#include "block_matrix/indexing.h"

#include "mp_tensors/compression.h"
#include "mp_tensors/contractions.h"

template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>
evolve(MPS<Matrix, SymmGroup> mps,
       block_matrix<Matrix, SymmGroup> const & op,
       std::size_t Mmax, double cutoff)
{
    std::size_t L = mps.length();
    
    for (int i = 0; i < 2; ++i) { // odd/even
        mps.canonize(i+1);
        block_matrix<Matrix, SymmGroup> v0, v1, t;
        
        for (std::size_t p = i; p < L-1; p += 2)
        {
//            cout << "Doing " << p << " " << p+1 << endl;
            
            mps[p].make_left_paired();
            mps[p+1].make_right_paired();
            
            gemm(mps[p].data(), mps[p+1].data(), v0);
            
            v1 = contraction::multiply_with_twosite(v0, op,
                                                    mps[p].row_dim(), mps[p+1].col_dim(),
                                                    mps[p].site_dim());
            
            compression::replace_two_sites(mps, Mmax, cutoff, v1, p);
            
            // move two to the right, if possible
            t = mps[p+1].normalize_left(SVD);
            if (p+2 < L) {
                mps[p+2].multiply_from_left(t);
                t = mps[p+2].normalize_left(SVD);
                if (p+3 < L) {
                    mps[p+3].multiply_from_left(t);
                }
            }
        }
        
//        cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << endl;
    }
    
    return mps;
}

#endif
