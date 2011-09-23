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
#include "utils/ambient_assert.h"
#include "utils/matrix_vector_traits.h"

#include "block_matrix/block_matrix.h"
#include "block_matrix/indexing.h"

#include <alps/numeric/real.hpp>
#include <alps/numeric/imag.hpp>

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
        C.insert_block(Matrix3(num_rows(A[k]), num_cols(B[matched_block])),
                       A.left_basis()[k].first, B.right_basis()[matched_block].first);
        gemm(A[k], B[matched_block], C[C.left_basis().position(A.left_basis()[k].first)]);
    }
}

template<class Tag> struct basis_eval;

template<> struct basis_eval<blas::NoTranspose>
{
    template<class Matrix, class SymmGroup> static Index<SymmGroup> const &
    first(block_matrix<Matrix, SymmGroup> const & m) { return m.left_basis(); }
    
    template<class Matrix, class SymmGroup> static Index<SymmGroup> const &
    second(block_matrix<Matrix, SymmGroup> const & m) { return m.right_basis(); }
};

template<> struct basis_eval<blas::Transpose>
{
    template<class Matrix, class SymmGroup> static Index<SymmGroup> const &
    first(block_matrix<Matrix, SymmGroup> const & m) { return m.right_basis(); }
    
    template<class Matrix, class SymmGroup> static Index<SymmGroup> const &
    second(block_matrix<Matrix, SymmGroup> const & m) { return m.left_basis(); }
};

template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup, class Tag1, class Tag2>
void gemm(block_matrix<Matrix1, SymmGroup> const & A, Tag1,
          block_matrix<Matrix2, SymmGroup> const & B, Tag2,
          block_matrix<Matrix3, SymmGroup> & C)
{
    C.clear();
    
    typedef typename SymmGroup::charge charge;
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        if (! basis_eval<Tag2>::first(B).has(basis_eval<Tag1>::second(A)[k].first))
            continue;
        
        std::size_t matched_block = basis_eval<Tag2>::first(B).position(basis_eval<Tag1>::second(A)[k].first);
        
        // avoid copying, use resize
        std::size_t s1 = result_size(A[k], Tag1(), B[matched_block], Tag2()).first;
        std::size_t s2 = result_size(A[k], Tag1(), B[matched_block], Tag2()).second;
        
        C.insert_block(Matrix3(s1, s2),
                       basis_eval<Tag1>::first(A)[k].first,
                       basis_eval<Tag2>::second(B)[matched_block].first);
        
        gemm(A[k], Tag1(), B[matched_block], Tag2(), C[C.left_basis().position(basis_eval<Tag1>::first(A)[k].first)]);
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
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(dynamic)
#endif
    for (std::size_t k = 0; k < loop_max; ++k)
        svd(M[k], U[k], V[k], S[k]);
    
    timer.end();
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void heev(block_matrix<Matrix, SymmGroup> const & M,
          block_matrix<Matrix, SymmGroup> & evecs,
          block_matrix<DiagMatrix, SymmGroup> & evals)
{
    static Timer timer("block_matrix EVD");
    timer.begin();

    evecs = block_matrix<Matrix, SymmGroup>(M.left_basis(), M.right_basis());
    evals = block_matrix<DiagMatrix, SymmGroup>(M.left_basis(), M.right_basis());
    std::size_t loop_max = M.n_blocks();
#ifndef MPI_PARALLEL
#pragma omp parallel for schedule(dynamic)
#endif
    for(std::size_t k = 0; k < loop_max; ++k)
        heev(M[k], evecs[k], evals[k]);

    timer.end();
}
//this function could return a void presently, check with Bela 

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd_truncate(block_matrix<Matrix, SymmGroup> const & M,
                  block_matrix<Matrix, SymmGroup> & U,
                  block_matrix<Matrix, SymmGroup> & V,
                  block_matrix<DiagMatrix, SymmGroup> & S,
                  double rel_tol, std::size_t Mmax,
                  bool verbose = true,
                  Logger * logger = NULL)
{   
  ambient_assert(false); 
  svd(M, U, V, S);
    
   //  Given the full SVD in each block (above), remove all singular values and corresponding rows/cols
   //  where the singular value is < rel_tol*max(S), where the maximum is taken over all blocks.
   //  Be careful to update the Index descriptions in the matrices to reflect the reduced block sizes
   //  (remove_rows/remove_columns methods for that)
    
    Index<SymmGroup> old_basis = S.left_basis();
    
    typename blas::associated_real_vector<Matrix>::type allS;
    
#ifdef MPI_PARALLEL
    ambient_assert(false);
#else
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        std::copy(S[k].elements().first, S[k].elements().second, std::back_inserter(allS));
#endif

#ifdef MPI_PARALLEL
  ambient_assert(false);
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
   ambient_assert(false);
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
    if (verbose /* && ! (old_basis == S.left_basis()) */) {
//        zout << "SVD performed a truncation: (cutoff = " << rel_tol << ")" << std::endl;
//        zout << old_basis << std::endl << S.left_basis() << std::endl;
        zout << "Sum: " << old_basis.sum_of_sizes() << " -> " << S.left_basis().sum_of_sizes() << std::endl;
//        zout << "Smallest SV kept: " << Scut / allS[0] << std::endl;
//        zout << "Truncated weight: " << truncated_weight << std::endl;
    }
    
//    if (logger != NULL) {
//        *logger << make_log("BondDimension", S.left_basis().sum_of_sizes());
//        *logger << make_log("TruncatedWeight", truncated_weight);
//        *logger << make_log("SmallestEV", Scut / allS[0]);
//    }
}

/*
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        blas::copy<Matrix>(S[k],allS);
    blas::sort<Matrix>(allS);  
   blas::reverse<Matrix>(allS); */

template<class Matrix, class DiagMatrix, class SymmGroup>
void heev_truncate(block_matrix<Matrix, SymmGroup> const & M,
                   block_matrix<Matrix, SymmGroup> & evecs,
                   block_matrix<DiagMatrix, SymmGroup> & evals,
                   double cutoff, std::size_t Mmax,
                   Logger & logger,
                   bool verbose = true)
{

    // very analogous to the above svd method
    assert( M.left_basis().sum_of_sizes() > 0 && M.right_basis().sum_of_sizes() > 0 );
    heev(M, evecs, evals);
    
    Index<SymmGroup> old_basis = evals.left_basis();
    
#ifdef MPI_PARALLEL
    size_t length = 0;
    size_t position = 0;
    for(std::size_t k = 0; k < evals.n_blocks(); ++k) 
        length += evals[k].size();
    std::vector<typename Matrix::value_type> allevals(length);
    for(std::size_t k = 0; k < evals.n_blocks(); ++k){
        blas::copy_after<Matrix>(allevals, position, evals[k]);
        position += evals[k].size();
    }
    ambient::playout(); // execution weight: 64
#else
    std::vector<double> allevals;
    for(std::size_t k = 0; k < evals.n_blocks(); ++k)
//        std::copy(evals[k].elements().first, evals[k].elements().second, std::back_inserter(allevals));
        for (typename DiagMatrix::const_element_iterator it = evals[k].elements().first; it != evals[k].elements().second; ++it) {
            assert( check_real(*it) );
            assert( alps::numeric::real(*it) > -1e-10 );
            allevals.push_back(alps::numeric::real(*it));
        }
#endif
    assert( allevals.size() > 0 );
    std::sort(allevals.begin(), allevals.end());
    std::reverse(allevals.begin(), allevals.end());
    /*std::cout << "All evals: ";

    for (std::size_t k = 0; k < allevals.size(); ++k) // LAUSANNE
        cout << allevals[k] << " ";
    cout << endl;*/

/* to do
*/
    
    double evalscut = cutoff * allevals[0];
    if (allevals.size() > Mmax)
        evalscut = std::max(evalscut, allevals[Mmax]);

    #ifndef MPI_PARALLEL
    double truncated_weight = std::accumulate(std::find_if(allevals.begin(), allevals.end(), boost::lambda::_1 < evalscut), allevals.end(), 0.0);
    truncated_weight /= std::accumulate(allevals.begin(), allevals.end(), 0.0);
    #endif

    #ifdef MPI_PARALLEL
    size_t* keeps = (size_t*)malloc(evals.n_blocks()*sizeof(size_t));
    for(size_t k = 0; k < evals.n_blocks(); ++k)
    {
        keeps[k] = num_rows(evals[k]);
        size_t* keep_ptr = keeps + k;
        ambient::push(ambient::associated_find_if_l, ambient::associated_find_if_c, evals[k].get_data(), evalscut, keep_ptr);
    }
    ambient::playout(); // execution weight 64
    #endif

    for (std::size_t k = 0; k < evals.n_blocks(); ++k)
    {
        #ifdef MPI_PARALLEL
        size_t keep = keeps[k];
        #else
        // FIXME: this is a stupid workaround
        std::vector<double> evals_k;
        for (typename DiagMatrix::const_element_iterator it = evals[k].elements().first; it != evals[k].elements().second; ++it)
            evals_k.push_back(alps::numeric::real(*it));
        size_t keep = std::find_if(evals_k.begin(), evals_k.end(),
                                   boost::lambda::_1 < evalscut)-evals_k.begin();
//        size_t keep = std::find_if(evals[k].elements().first, evals[k].elements().second,
//                                        boost::lambda::_1 < evalscut)-evals[k].elements().first;
        #endif

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
   
    if (verbose /* && ! (old_basis == evals.left_basis())*/ ) {
//        zout << "heev_truncate performed: (cutoff = " << cutoff << ")" << std::endl;
//        zout << old_basis << std::endl << evals.left_basis() << std::endl;
        zout << "Sum: " << old_basis.sum_of_sizes() << " -> " << evals.left_basis().sum_of_sizes() << std::endl;
//        zout << "Smallest EV kept: " << evalscut / allevals[0] << std::endl;
    }
    
    logger << make_log("BondDimension", evals.left_basis().sum_of_sizes());
    #ifndef MPI_PARALLEL
    logger << make_log("TruncatedWeight", truncated_weight);
    #endif
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
        m[k] = sqrt(m[k]);

    return m;
}

// Is it really useful?
template <class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> exp (block_matrix<Matrix, SymmGroup> const & M, typename Matrix::value_type const & alpha = 1)
{
    block_matrix<Matrix, SymmGroup> N, tmp, res;
    block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
    
    heev(M, N, S);
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        S[k] = exp(alpha*S[k]);
    gemm(N, S, tmp);
    gemm(tmp, conjugate_transpose(N), res);
    
    return res;
}

template <class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> op_exp (Index<SymmGroup> const & phys,
                                        block_matrix<Matrix, SymmGroup> M,
                                        typename Matrix::value_type const & alpha = 1)
{
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            M(it_c->first, it_c->first) = exp(M(it_c->first, it_c->first), alpha);
        else
            M.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return M;
}

template<class Matrix, class SymmGroup>
void op_kron(Index<SymmGroup> const & phys,
             block_matrix<Matrix, SymmGroup> const & A,
             block_matrix<Matrix, SymmGroup> const & B,
             block_matrix<Matrix, SymmGroup> & C)
{
    C = block_matrix<Matrix, SymmGroup>();
    
    Index<SymmGroup> const & left_basis = phys;
    Index<SymmGroup> const & right_basis = phys;
    
    ProductBasis<SymmGroup> pb_left(left_basis, left_basis);
    ProductBasis<SymmGroup> pb_right(right_basis, right_basis);
    
    for (int i=0; i<A.n_blocks(); ++i) {
        for (int j=0; j<B.n_blocks(); ++j) {
            typename SymmGroup::charge new_right = SymmGroup::fuse(A.right_basis()[i].first, B.right_basis()[j].first);
            typename SymmGroup::charge new_left = SymmGroup::fuse(A.left_basis()[i].first, B.left_basis()[j].first);
            
            
            Matrix tmp(pb_left.size(A.left_basis()[i].first, B.left_basis()[j].first),
                       pb_right.size(A.right_basis()[i].first, B.right_basis()[j].first),
                       0);

            for (int l1 = 0; l1 < A.left_basis()[i].second; ++l1)
                for (int l2 = 0; l2 < B.left_basis()[j].second; ++l2)
                    for (int r1 = 0; r1 < A.right_basis()[i].second; ++r1)
                        for (int r2 = 0; r2 < B.right_basis()[j].second; ++r2)
                            tmp(pb_left(A.left_basis()[i].first,
                                        B.left_basis()[j].first)+l1*B.left_basis()[j].second+l2,
                                pb_right(A.right_basis()[i].first,
                                         B.right_basis()[j].first)+r1*B.right_basis()[j].second+r2) = A[i](l1, r1) * B[j](l2, r2);
            
            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
    
}

#endif
