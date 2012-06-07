/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef EVOLVE_H
#define EVOLVE_H

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
       std::size_t Mmax, double cutoff, Logger * logger = NULL)
{
    std::size_t L = mps.length();
    
    for (int i = 0; i < 2; ++i) { // odd/even
        mps.canonize(i+1);
        block_matrix<Matrix, SymmGroup> v0, v1, t;
        
        for (std::size_t p = i; p < L-1; p += 2)
        {
//            maquis::cout << "Doing " << p << " " << p+1 << std::endl;
            
            mps[p].make_left_paired();
            mps[p+1].make_right_paired();
            
            gemm(mps[p].data(), mps[p+1].data(), v0);
            
            v1 = contraction::multiply_with_twosite(v0, op,
                                                    mps[p].row_dim(), mps[p+1].col_dim(),
                                                    mps[p].site_dim());
            
            compression::replace_two_sites_l2r(mps, Mmax, cutoff, v1, p, logger);
            
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
        
//        maquis::cout << "Norm loss " << i << ": " << trace(t)-1.0 << " " << -log(trace(t)) << std::endl;
    }
    
    return mps;
}

// Same function, but working with different matrix on each bond.
// map_op should already contain non overlapping terms.
template<class Matrix, class SymmGroup>
void
evolve_l2r(MPS<Matrix, SymmGroup> & mps,
           std::vector<block_matrix<Matrix, SymmGroup> > const & ops,
           std::vector<long> const & idx,
           int pfirst, std::size_t Mmax, double cutoff, Logger * logger = NULL)
{
    assert(mps.length() == idx.size());
    std::size_t L = mps.length();
    
    MPS<Matrix, SymmGroup> const& constmps = mps;
    
    if (mps.canonization() != pfirst + 1)
        mps.canonize(pfirst + 1);

    for (std::size_t p = pfirst; p <= L-1; p += 2)
    {
        if (idx[p] != -1)
        {
            constmps[p].make_left_paired();
            constmps[p+1].make_right_paired();
            
            block_matrix<Matrix, SymmGroup> v0, v1;
            gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
            
            v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                    constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                    constmps[p].site_dim());
            compression::replace_two_sites_l2r(mps, Mmax, cutoff, v1, p, logger);
        }
        mps.move_normalization_l2r(p+1, p+3);
    }
    mps.canonization(true);
    assert(mps.canonization() == 0);
    // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
}

template<class Matrix, class SymmGroup>
void
evolve_r2l(MPS<Matrix, SymmGroup> & mps,
           std::vector<block_matrix<Matrix, SymmGroup> > const & ops,
           std::vector<long> const & idx,
           int pfirst, std::size_t Mmax, double cutoff, Logger * logger = NULL)
{
    assert(mps.length() == idx.size());
    std::size_t L = mps.length();

    MPS<Matrix, SymmGroup> const& constmps = mps;
    
    int startpos = std::min(L-2-(L-pfirst)%2, L-2);
    if (mps.canonization() != startpos)
        mps.canonize(startpos);

    for (int p = std::min(L-2-(L-pfirst)%2, L-2); p >= pfirst; p -= 2)
    {
        if (idx[p] != -1)
        {
            mps[p].make_left_paired();
            mps[p+1].make_right_paired();

            block_matrix<Matrix, SymmGroup> v0, v1;
            gemm(constmps[p].data(), constmps[p+1].data(), v0); // outer product of two sites
            
            v1 = contraction::multiply_with_twosite(v0, ops[idx[p]],
                                                    constmps[p].row_dim(), constmps[p+1].col_dim(),
                                                    constmps[p].site_dim());
            compression::replace_two_sites_r2l(mps, Mmax, cutoff, v1, p, logger);
        }
        mps.move_normalization_r2l(p, std::max(static_cast<long>(p)-2,0L));
    }
    
    
    mps.canonization(true);
    assert(mps.canonization() == mps.length()-1);
    // maquis::cout << "Norm loss " << i << ": " << trace(t) << " " << -log(trace(t)) << std::endl;
}

#endif
