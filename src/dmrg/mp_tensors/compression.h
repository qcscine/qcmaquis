/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef COMPRESSION_H
#define COMPRESSION_H

#include "utils/zout.hpp"

#include "mp_tensors/mpstensor.h"
#include "mp_tensors/mpotensor.h"

#include "mp_tensors/reshapes.h"
#include "block_matrix/indexing.h"

struct compression {
    template<class Matrix, class SymmGroup>
    void static
    replace_two_sites(MPS<Matrix, SymmGroup> & mps,
                      std::size_t Mmax, double cutoff,
                      block_matrix<Matrix, SymmGroup> const & t,
                      std::size_t p,
                      Logger * logger = NULL)
    {
        block_matrix<Matrix, SymmGroup> u, v;
        
        typedef typename blas::associated_real_diagonal_matrix<Matrix>::type dmt;

        block_matrix<dmt, SymmGroup> s;
        
        svd_truncate(t, u, v, s,
                     cutoff, Mmax, false, logger);
        
        mps[p].make_left_paired();
        mps[p+1].make_right_paired();
        
        mps[p].data() = u;
        mps[p].right_i = u.right_basis();
        
        gemm(s, v, u);
        mps[p+1].data() = u;
        mps[p+1].left_i = u.left_basis();
    }
    
    template<class Matrix, class SymmGroup>
    static void compress_two_sites(MPS<Matrix, SymmGroup> & mps,
                                   std::size_t Mmax, double cutoff,
                                   std::size_t p,
                                   Logger * logger = NULL)
    {
        block_matrix<Matrix, SymmGroup> t;
        
        mps[p].make_left_paired();
        mps[p+1].make_right_paired();
        
        gemm(mps[p].data(), mps[p+1].data(), t);
        
        replace_two_sites(mps, Mmax, cutoff, t, p, logger);
    }
        
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static l2r_compress(MPS<Matrix, SymmGroup> mps,
                        std::size_t Mmax, double cutoff,
                        Logger * logger = NULL,
                        bool verbose = false)
    {
        std::size_t L = mps.length();
        std::vector<double> ret;
        
        block_matrix<Matrix, SymmGroup> t;
        
        mps.canonize(1);
        
        if (verbose) cout << "Compressing @ ";
        for (std::size_t p = 1; p < L; ++p)
        {
            if (verbose) {
                cout << p << " ";
                cout.flush();
            }
            
            compress_two_sites(mps, Mmax, cutoff, p-1, logger);
            
            t = mps[p].normalize_left(SVD);
            
            if (p+1 < L)
                mps[p+1].multiply_from_left(t);
            else
                zout << "Norm reduction: " << trace(t) << endl;
        }
        
        return mps;
    }
    
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup>
    static r2l_compress(MPS<Matrix, SymmGroup> mps,
                        std::size_t Mmax, double cutoff,
                        Logger * logger = NULL,
                        bool verbose = false)
    {
        std::size_t L = mps.length();
        std::vector<double> ret;
        
        block_matrix<Matrix, SymmGroup> t;
        
        mps.canonize(L-1);
        
        if (verbose) cout << "Compressing @ ";
        for (std::size_t p = L-1; p > 0; --p)
        {
            if (verbose) {
                cout << p << " ";
                cout.flush();
            }
            
            compress_two_sites(mps, Mmax, cutoff, p-1, logger);
            
            t = mps[p-1].normalize_right(SVD);
            
            if (p > 1)
                mps[p-2].multiply_from_right(t);
            else
                zout << "Norm reduction: " << trace(t) << endl;
        }
        
        return mps;
    }
};

#endif
