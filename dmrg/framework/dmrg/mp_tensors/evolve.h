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

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/indexing.h"

#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/contractions.h"

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
        
//        cout << "Norm loss " << i << ": " << trace(t)-1.0 << " " << -log(trace(t)) << endl;
    }
    
    return mps;
}

// Same function, but working with different matrix on each bond.
// map_op should already contain non overlapping terms.
template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup>
evolve(MPS<Matrix, SymmGroup> mps,
       std::map<std::size_t, block_matrix<Matrix, SymmGroup> > const & map_op,
       std::size_t Mmax, double cutoff)
{
    std::size_t L = mps.length();
    
    int i = (map_op.begin()->first % 2); // odd/even
    
    mps.canonize(i+1);
    block_matrix<Matrix, SymmGroup> v0, v1, t;
    
    for (std::size_t p = i; p < L-1; p += 2)
    {
        
        mps[p].make_left_paired();
        mps[p+1].make_right_paired();
        
        if (map_op.count(p) > 0) {
            gemm(mps[p].data(), mps[p+1].data(), v0);
            
            v1 = contraction::multiply_with_twosite(v0, map_op.find(p)->second,
                                                    mps[p].row_dim(), mps[p+1].col_dim(),
                                                    mps[p].site_dim());
            
            compression::replace_two_sites(mps, Mmax, cutoff, v1, p);
            
        }
        
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
    
    return mps;
}

#endif
