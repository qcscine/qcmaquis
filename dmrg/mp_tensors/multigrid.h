/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MULTIGRID_H
#define MULTIGRID_H

#include "utils/zout.hpp"

#include "mp_tensors/mpstensor.h"
#include "block_matrix/indexing.h"

struct multigrid {
    
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static restriction ();
    
    template<class Matrix, class SymmGroup>
    void
    static extension (MPS<Matrix, SymmGroup> mps_small,
                      MPS<Matrix, SymmGroup> & mps_large,
                       std::size_t Mmax, double cutoff)
    {
        std::size_t L1 = mps_small.length();
        std::size_t L2 = mps_large.length();
        assert(L2 == 2*L1);
                        
        for (std::size_t p = 1; p < L1; ++p)
        {
            Index<SymmGroup> phys = mps_small[p].site_dim();
            std::vector<block_matrix<Matrix, SymmGroup> > mlist = mps_small[p].to_list();
            
            std::vector<block_matrix<Matrix, SymmGroup> > left(mlist.size()), right(mlist.size());
            for (int i=0; i<mlist.size(); ++i)
            {
                block_matrix<Matrix, SymmGroup> U, V;
                block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Ssqrt;
                svd(mlist[i], U, V, S);
                Ssqrt = sqrt(S);
                gemm(U, Ssqrt, left[i].second);
                gemm(Ssqrt, V, right[i].second);
            }
            
            mps_large[2*p] = MPSTensor(phys, left);
            mps_large[2*p+1] = MPSTensor(phys, right);
        }
        
    }
    
};

#endif
