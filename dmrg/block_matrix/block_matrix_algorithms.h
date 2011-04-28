/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BLOCK_MATRIX_ALGORITHMS_H
#define BLOCK_MATRIX_ALGORITHMS_H

#include "utils/zout.hpp"
#include "utils/logger.h"
#include "utils/timings.h"

#include "block_matrix/block_matrix.h"
// some example functions
template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm(block_matrix<Matrix1, SymmGroup> const & A,
          block_matrix<Matrix2, SymmGroup> const & B,
          block_matrix<Matrix3, SymmGroup> & C)
{
    C.clear();
    
    typedef typename SymmGroup::charge charge;
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        if (! B.left_basis().has(A.right_basis()[k].first))
            continue;
        
        std::size_t matched_block = B.left_basis().position(A.right_basis()[k].first);
        
        // avoid copying, use resize
        C.insert_block(Matrix3(),
                       A.left_basis()[k].first, B.right_basis()[matched_block].first);
        C.resize_block(A.left_basis()[k].first, B.right_basis()[matched_block].first,
                       num_rows(A[k]), num_cols(B[matched_block]));
        gemm(A[k], B[matched_block], C[C.left_basis().position(A.left_basis()[k].first)]);
    }
}

template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void pgemm(block_matrix<Matrix1, SymmGroup> const & A,
          block_matrix<Matrix2, SymmGroup> const & B,
          block_matrix<Matrix3, SymmGroup> & C)
{
    C.clear();
    
    typedef typename SymmGroup::charge charge;
    std::size_t loop_max = A.n_blocks();
    for (std::size_t k = 0; k < loop_max; ++k) {
        if (! B.left_basis().has(A.right_basis()[k].first))
            continue;

        std::size_t matched_block = B.left_basis().position(A.right_basis()[k].first);

        // avoid copying, use resize
        C.insert_block(Matrix3(),
                       A.left_basis()[k].first, B.right_basis()[matched_block].first);
        C.resize_block(A.left_basis()[k].first, B.right_basis()[matched_block].first,
                       num_rows(A[k]), num_cols(B[matched_block]));
    }

#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < loop_max; ++k) {
        if (! B.left_basis().has(A.right_basis()[k].first))
            continue;
        
        std::size_t matched_block = B.left_basis().position(A.right_basis()[k].first);
        
        // avoid copying, use resize
        gemm(A[k], B[matched_block], C[C.left_basis().position(A.left_basis()[k].first)]);
    }
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd(block_matrix<Matrix, SymmGroup> const & M,
         block_matrix<Matrix, SymmGroup> & U,
         block_matrix<Matrix, SymmGroup> & V,
         block_matrix<DiagMatrix, SymmGroup> & S)
{
    static Timer timer("block_matrix SVD");
    timer.begin();
    
    Index<SymmGroup> r = M.left_basis(), c = M.right_basis(), m = M.left_basis();
    for (std::size_t i = 0; i < M.n_blocks(); ++i)
        m[i].second = std::min(r[i].second, c[i].second);
    
    U = block_matrix<Matrix, SymmGroup>(r, m);
    V = block_matrix<Matrix, SymmGroup>(m, c);
    S = block_matrix<DiagMatrix, SymmGroup>(m, m);
    
    std::size_t loop_max = M.n_blocks();
#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < loop_max; ++k)
        svd(M[k], U[k], V[k], S[k]);
    
    timer.end();
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void syev(block_matrix<Matrix, SymmGroup> const & M,
          block_matrix<Matrix, SymmGroup> & evecs,
          block_matrix<DiagMatrix, SymmGroup> & evals)
{
    static Timer timer("block_matrix EVD");
    timer.begin();

    evecs = block_matrix<Matrix, SymmGroup>(M.left_basis(), M.right_basis());
    evals = block_matrix<DiagMatrix, SymmGroup>(M.left_basis(), M.right_basis());
    std::size_t loop_max = M.n_blocks();
#pragma omp parallel for schedule(dynamic)
    for (std::size_t k = 0; k < loop_max; ++k)
        syev(M[k], evecs[k], evals[k]);

    timer.end();
}
//this function could return a void presently, check with Bela 

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd_truncate(block_matrix<Matrix, SymmGroup> const & M,
                  block_matrix<Matrix, SymmGroup> & U,
                  block_matrix<Matrix, SymmGroup> & V,
                  block_matrix<DiagMatrix, SymmGroup> & S,
                  double rel_tol, std::size_t Mmax,
                  bool verbose = true)
{   
  assert(false); 
  svd(M, U, V, S);
    
   //  Given the full SVD in each block (above), remove all singular values and corresponding rows/cols
   //  where the singular value is < rel_tol*max(S), where the maximum is taken over all blocks.
   //  Be careful to update the Index descriptions in the matrices to reflect the reduced block sizes
   //  (remove_rows/remove_columns methods for that)
    
    Index<SymmGroup> old_basis = S.left_basis();
    
    
    typename blas::associated_vector<Matrix>::type allS;
#ifdef MPI_PARALLEL
 assert(false);
#else

    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        std::copy(S[k].elements().first, S[k].elements().second, std::back_inserter(allS));
#endif


#ifdef MPI_PARALLEL
 assert(false);
  double Scut;
#else
    std::sort(allS.begin(), allS.end());
    std::reverse(allS.begin(), allS.end());
    double Scut = rel_tol * allS[0];
    if (allS.size() > Mmax)
        Scut = std::max(Scut, allS[Mmax]);
    double truncated_weight = std::accumulate(std::find_if(allS.begin(), allS.end(), boost::lambda::_1 < Scut), allS.end(), 0.0);
    truncated_weight /= std::accumulate(allS.begin(), allS.end(), 0.0);
   #endif 
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
    {
#ifdef MPI_PARALLEL
   assert(false);
   int keep = 0;
#else
        std::size_t keep = std::find_if(S[k].elements().first, S[k].elements().second,
                                        boost::lambda::_1 < Scut)-S[k].elements().first;
#endif 
       if (keep >= num_rows(S[k]))
            continue;
        
    //  hack for now 
//        keep = std::max(keep, std::size_t(1));
        
        if (keep == 0) {
            S.remove_block(S.left_basis()[k].first,
                           S.right_basis()[k].first);
            U.remove_block(U.left_basis()[k].first,
                           U.right_basis()[k].first);
            V.remove_block(V.left_basis()[k].first,
                           V.right_basis()[k].first);
            --k; // everything gets shifted, to we have to look into the same k again
        } else {
            S.resize_block(S.left_basis()[k].first,
                           S.right_basis()[k].first,
                           keep, keep);
            
            U.resize_block(U.left_basis()[k].first,
                           U.right_basis()[k].first,
                           U.left_basis()[k].second,
                           keep);
            
            V.resize_block(V.left_basis()[k].first,
                           V.right_basis()[k].first,
                           keep,
                           V.right_basis()[k].second);
        }
    }
   /* 
    if (verbose && ! (old_basis == S.left_basis()) ) {
        zout << "SVD performed a truncation: (cutoff = " << rel_tol << ")" << std::endl;
        zout << old_basis << std::endl << S.left_basis() << std::endl;
        zout << "Sum: " << old_basis.sum_of_sizes() << " -> " << S.left_basis().sum_of_sizes() << std::endl;
        zout << "Smallest SV kept: " << Scut / allS[0] << std::endl;
        zout << "Truncated weight: " << truncated_weight << std::endl;
    }
    */
}

/*
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        blas::copy<Matrix>(S[k],allS);
    blas::sort<Matrix>(allS);  
   blas::reverse<Matrix>(allS); */

template<class Matrix, class DiagMatrix, class SymmGroup>
void syev_truncate(block_matrix<Matrix, SymmGroup> const & M,
                                             block_matrix<Matrix, SymmGroup> & evecs,
                                             block_matrix<DiagMatrix, SymmGroup> & evals,
                                             double cutoff, std::size_t Mmax,
                                             Logger & logger,
                                             bool verbose = true)
{

 // assert(false); 
    // very analogous to the above svd method
    syev(M, evecs, evals);
    
    Index<SymmGroup> old_basis = evals.left_basis();
    
    //std::vector<double> allevals;
    typename blas::associated_vector<Matrix>::type allevals;
assert(false);
/* to do
     for (std::size_t k = 0; k < evals.n_blocks(); ++k)
        std::copy(evals[k].elements().first, evals[k].elements().second, std::back_inserter(allevals));
*/
    std::sort(allevals.begin(), allevals.end());
    std::reverse(allevals.begin(), allevals.end());

    
    double evalscut = cutoff * allevals[0];
    if (allevals.size() > Mmax)
        evalscut = std::max(evalscut, allevals[Mmax]);
    double truncated_weight = std::accumulate(std::find_if(allevals.begin(), allevals.end(), boost::lambda::_1 < evalscut), allevals.end(), 0.0);
    truncated_weight /= std::accumulate(allevals.begin(), allevals.end(), 0.0);
    
    for (std::size_t k = 0; k < evals.n_blocks(); ++k)
  {
        int keep = std::find_if(evals[k].elements().first, evals[k].elements().second,
                                        boost::lambda::_1 < evalscut)-evals[k].elements().first;


        if (keep >= num_rows(evals[k]))
            continue;
        
        // uncomment if you want to keep all blocks 
//        keep = std::max(keep, std::size_t(1));
        
        if (keep == 0) {
            evals.remove_block(evals.left_basis()[k].first,
                               evals.right_basis()[k].first);
            evecs.remove_block(evecs.left_basis()[k].first,
                               evecs.right_basis()[k].first);
            --k; // everything gets shifted, to we have to look into the same k again
        } else {
            evals.resize_block(evals.left_basis()[k].first,
                               evals.right_basis()[k].first,
                               keep, keep);
            
            evecs.resize_block(evecs.left_basis()[k].first,
                               evecs.right_basis()[k].first,
                               evecs.left_basis()[k].second,
                               keep);
        }
    }
   /* 
    if (verbose && ! (old_basis == evals.left_basis()) ) {
        zout << "syev_truncate performed: (cutoff = " << cutoff << ")" << std::endl;
        zout << old_basis << std::endl << evals.left_basis() << std::endl;
        zout << "Sum: " << old_basis.sum_of_sizes() << " -> " << evals.left_basis().sum_of_sizes() << std::endl;
        zout << "Smallest EV kept: " << evalscut / allevals[0] << std::endl;
    }
    */
    logger << make_log("BondDimension", evals.left_basis().sum_of_sizes());
    logger << make_log("TruncatedWeight", truncated_weight);
    logger << make_log("SmallestEV", evalscut / allevals[0]);
}



template<class Matrix, class SymmGroup>
void qr(block_matrix<Matrix, SymmGroup> & M,
        block_matrix<Matrix, SymmGroup> & Q,
        block_matrix<Matrix, SymmGroup> & R)
{
    /* thin QR in each block */
    
    Index<SymmGroup> m = M.left_basis(), n = M.right_basis();
    
    Q = block_matrix<Matrix, SymmGroup>(m,n);
    R = block_matrix<Matrix, SymmGroup>(n,n);
    
    for (std::size_t k = 0; k < M.n_blocks(); ++k)
        qr(M[k], Q[k], R[k]);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> transpose(block_matrix<Matrix, SymmGroup> m)
{
    m.inplace_transpose();
    return m;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> conjugate(block_matrix<Matrix, SymmGroup> m)
{
    m.inplace_conjugate();
    return m;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> conjugate_transpose(block_matrix<Matrix, SymmGroup> m)
{
    m.inplace_transpose();
    m.inplace_conjugate();
    return m;
}

template<class Matrix, class SymmGroup>
typename Matrix::value_type trace(block_matrix<Matrix, SymmGroup> const & m)
{
    return m.trace();
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> adjoin(block_matrix<Matrix, SymmGroup> const & m)
{
    block_matrix<Matrix, SymmGroup> ret;
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        ret.insert_block(m[k],
                         -m.left_basis()[k].first,
                         -m.right_basis()[k].first);
    return ret;
}

template<class Matrix, class SymmGroup, class Generator>
void generate(block_matrix<Matrix, SymmGroup> & m, Generator & g)
{
    m.generate(g);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> identity_matrix(Index<SymmGroup> const & size)
{
    block_matrix<Matrix, SymmGroup> ret(size, size);
    for (std::size_t k = 0; k < ret.n_blocks(); ++k)
        ret[k] = Matrix::identity_matrix(size[k].second);
    return ret;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> sqrt(block_matrix<Matrix, SymmGroup>  m)
{
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
    {
#ifdef MPI_PARALLEL
        sqrt(m[k]);
#else
        m[k] = sqrt(m[k]);
#endif
    }

    return m;
}

#endif
