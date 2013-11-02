/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef BLOCK_MATRIX_ALGORITHMS_H
#define BLOCK_MATRIX_ALGORITHMS_H

#include "dmrg/utils/logger.h"
#include "dmrg/utils/utils.hpp"
#include "utils/timings.h"
#include "utils/traits.hpp"
#include "utils/bindings.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"

#include <boost/lambda/lambda.hpp>
#include <boost/function.hpp>

#include "dmrg/utils/parallel_for.hpp"

struct truncation_results {
    std::size_t bond_dimension;     // new bond dimension
    double      truncated_weight;   // sum of discarded eigenvalues (square of singuler values)
    double      truncated_fraction; // sum of discarded singular values
    double      smallest_ev;        // smallest eigenvalue kept

    truncation_results() { }
    
    truncation_results(std::size_t m, double tw, double tf, double se)
    : bond_dimension(m)
    , truncated_weight(tw)
    , truncated_fraction(tf)
    , smallest_ev(se)
    { }
};

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
        
        std::size_t new_block = C.insert_block(new Matrix3(num_rows(A[k]), num_cols(B[matched_block])),
                                               A.left_basis()[k].first, B.right_basis()[matched_block].first);
        gemm(A[k], B[matched_block], C[new_block]);
    }
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd(block_matrix<Matrix, SymmGroup> const & M,
         block_matrix<Matrix, SymmGroup> & U,
         block_matrix<Matrix, SymmGroup> & V,
         block_matrix<DiagMatrix, SymmGroup> & S)
{
    Index<SymmGroup> r = M.left_basis(), c = M.right_basis(), m = M.left_basis();
    for (std::size_t i = 0; i < M.n_blocks(); ++i)
        m[i].second = std::min(r[i].second, c[i].second);
    
    U = block_matrix<Matrix, SymmGroup>(r, m);
    V = block_matrix<Matrix, SymmGroup>(m, c);
    S = block_matrix<DiagMatrix, SymmGroup>(m, m);
    std::size_t loop_max = M.n_blocks();
    
    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        svd(M[k], U[k], V[k], S[k]);
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void heev(block_matrix<Matrix, SymmGroup> const & M,
          block_matrix<Matrix, SymmGroup> & evecs,
          block_matrix<DiagMatrix, SymmGroup> & evals)
{

    evecs = block_matrix<Matrix, SymmGroup>(M.left_basis(), M.right_basis());
    evals = block_matrix<DiagMatrix, SymmGroup>(M.left_basis(), M.right_basis());
    std::size_t loop_max = M.n_blocks();

    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        heev(M[k], evecs[k], evals[k]);
}
    
#ifdef AMBIENT
template<class Matrix, class DiagMatrix, class SymmGroup>
void svd_merged(block_matrix<Matrix, SymmGroup> const & M,
                block_matrix<Matrix, SymmGroup> & U,
                block_matrix<Matrix, SymmGroup> & V,
                block_matrix<DiagMatrix, SymmGroup> & S)
{
    Index<SymmGroup> r = M.left_basis(), c = M.right_basis(), m = M.left_basis();
    for (std::size_t i = 0; i < M.n_blocks(); ++i)
        m[i].second = std::min(r[i].second, c[i].second);
    
    U = block_matrix<Matrix, SymmGroup>(r, m);
    V = block_matrix<Matrix, SymmGroup>(m, c);
    S = block_matrix<DiagMatrix, SymmGroup>(m, m);
    std::size_t loop_max = M.n_blocks();

    #ifdef AMBIENT_REBALANCED_SVD
    static ambient::synctime timer("SVD ONLY TIME\n");
    static ambient::synctime timert("SVD TRANSFER/MERGE TIME\n");
    timert.begin();

    // calculating complexities of the svd calls
    size_t np = ambient::channel.wk_dim();
    double total = 0;
    std::vector< std::pair<double,size_t> > complexities(loop_max);
    for(size_t i = 0; i < loop_max; ++i){
        double rows = num_rows(M[i]);
        double cols = num_cols(M[i]);
        complexities[i] = std::make_pair(2*rows*rows*cols + 4*cols*cols*cols, i);
        total += complexities[i].first;
    }
    total /= np;

    std::sort(complexities.begin(), complexities.end(), [](const std::pair<double,size_t>& a, const std::pair<double,size_t>& b){ return a.first < b.first; });
    std::vector<std::pair<double,size_t> > workloads(np, std::make_pair(0,0));

    // filling the workload with smallest local matrices first
    for(size_t p = 0; p < np; ++p){
        workloads[p].second = p;
        locale l(p);
        for(size_t i = 0; i < loop_max; ++i){
            size_t k = complexities[i].second;
            if(M[k][0].core->current->owner != p && (p != ambient::rank() || M[k][0].core->current->owner != -1)) continue;
            if(workloads[p].first + complexities[i].first >= total) break;
            merge(M[k]); 
            workloads[p].first += complexities[i].first; 
            complexities[i].first = 0;
        }
    }

    // rebalancing using difference with average
    for(size_t p = 0; p < np; ++p){
        locale l(p);
        for(size_t i = 0; i < loop_max; ++i){
            size_t k = complexities[i].second;
            if(complexities[i].first == 0) continue;
            if(workloads[p].first + complexities[i].first >= total) break;
            merge(M[k]); 
            workloads[p].first += complexities[i].first; 
            complexities[i].first = 0;
        }
    }

    std::sort(workloads.begin(), workloads.end(), [](const std::pair<double,size_t>& a, const std::pair<double,size_t>& b){ return a.first < b.first; });

    // placing the rest according to sorted workload
    size_t p = 0;
    for(int i = loop_max-1; i >= 0; --i){
        if(complexities[i].first == 0) continue;
        size_t k = complexities[i].second;
        int owner = workloads[p++].second;
        locale l(owner);
        merge(M[k]); 
        p %= np;
    }
    timert.end();

    #ifdef AMBIENT_PARALLEL_MKL
    ambient::mkl_set_num_threads(AMBIENT_MKL_NUM_THREADS);
    #endif
    timer.begin();
    for(size_t k = 0; k < loop_max; ++k){
        int owner = M[k][0].core->current->owner;
        if(owner == -1) owner = ambient::rank();
        locale l(owner);
        svd_merged(M[k], U[k], V[k], S[k]);
    }
    timer.end();
    #ifdef AMBIENT_PARALLEL_MKL
    ambient::mkl_set_num_threads(1);
    #endif
    #else
    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        svd_merged(M[k], U[k], V[k], S[k]);
    #endif
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void heev_merged(block_matrix<Matrix, SymmGroup> const & M,
                 block_matrix<Matrix, SymmGroup> & evecs,
                 block_matrix<DiagMatrix, SymmGroup> & evals)
{

    evecs = block_matrix<Matrix, SymmGroup>(M.left_basis(), M.right_basis());
    evals = block_matrix<DiagMatrix, SymmGroup>(M.left_basis(), M.right_basis());
    std::size_t loop_max = M.n_blocks();

    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        heev_merged(M[k], evecs[k], evals[k]);
}
#endif

template <class T>
typename maquis::traits::real_type<T>::type gather_real_pred(T const & val)
{
    assert( check_real(val) );
    assert( maquis::real(val) > -1e-10 );
    return maquis::real(val);
}

template<class DiagMatrix, class SymmGroup>
void estimate_truncation(block_matrix<DiagMatrix, SymmGroup> const & evals, 
                         size_t Mmax, double cutoff, size_t* keeps, 
                         double & truncated_fraction, double & truncated_weight, double & smallest_ev)
{ // to be parallelized later (30.04.2012)
    typedef typename DiagMatrix::value_type value_type;

    size_t length = 0;
    for(std::size_t k = 0; k < evals.n_blocks(); ++k){
        length += num_rows(evals[k]);
    }
    
    typedef std::vector<typename maquis::traits::real_type<value_type>::type > real_vector_t;
    real_vector_t allevals(length);
#ifdef AMBIENT
    locale_shared i;
    for(std::size_t k = 0; k < evals.n_blocks(); ++k){
        ambient::numeric::migrate(const_cast<DiagMatrix&>(evals[k])[0]);
    }
    ambient::sync();
#endif
    
    std::size_t position = 0;
    for(std::size_t k = 0; k < evals.n_blocks(); ++k){
        std::transform(evals[k].diagonal().first, evals[k].diagonal().second, allevals.begin()+position, gather_real_pred<value_type>);
        position += num_rows(evals[k]);
    }
    
    assert( allevals.size() > 0 );
    std::sort(allevals.begin(), allevals.end());
    std::reverse(allevals.begin(), allevals.end());
    
    double evalscut = cutoff * allevals[0];
    
    if (allevals.size() > Mmax)
        evalscut = std::max(evalscut, allevals[Mmax]);
    smallest_ev = evalscut / allevals[0];
    
    truncated_fraction = 0.0; truncated_weight = 0.0;
    for (typename real_vector_t::const_iterator it = std::find_if(allevals.begin(), allevals.end(), boost::lambda::_1 < evalscut);
         it != allevals.end(); ++it) {
        truncated_fraction += *it;        
        truncated_weight += (*it)*(*it);
    }
    truncated_fraction /= std::accumulate(allevals.begin(), allevals.end(), 0.0);
    truncated_weight /= std::accumulate(allevals.begin(), allevals.end(), 0.0,  boost::lambda::_1 + boost::lambda::_2 *boost::lambda::_2);
    
    for(std::size_t k = 0; k < evals.n_blocks(); ++k){
        real_vector_t evals_k(num_rows(evals[k]));
        std::transform(evals[k].diagonal().first, evals[k].diagonal().second, evals_k.begin(), gather_real_pred<value_type>);
        keeps[k] = std::find_if(evals_k.begin(), evals_k.end(), boost::lambda::_1 < evalscut)-evals_k.begin();
    }
}


template<class Matrix, class DiagMatrix, class SymmGroup>
truncation_results svd_truncate(block_matrix<Matrix, SymmGroup> const & M,
                                block_matrix<Matrix, SymmGroup> & U,
                                block_matrix<Matrix, SymmGroup> & V,
                                block_matrix<DiagMatrix, SymmGroup> & S,
                                double rel_tol, std::size_t Mmax,
                                bool verbose = true)
{ 
    assert( M.left_basis().sum_of_sizes() > 0 && M.right_basis().sum_of_sizes() > 0 );
    #ifdef AMBIENT
    svd_merged(M, U, V, S);
    #else
    svd(M, U, V, S);
    #endif
    
    Index<SymmGroup> old_basis = S.left_basis();
    size_t* keeps = new size_t[S.n_blocks()];
    double truncated_fraction, truncated_weight, smallest_ev;
    //  Given the full SVD in each block (above), remove all singular values and corresponding rows/cols
    //  where the singular value is < rel_tol*max(S), where the maximum is taken over all blocks.
    //  Be careful to update the Index descriptions in the matrices to reflect the reduced block sizes
    //  (remove_rows/remove_cols methods for that)
    estimate_truncation(S, Mmax, rel_tol, keeps, truncated_fraction, truncated_weight, smallest_ev);
     
    for ( int k = S.n_blocks() - 1; k >= 0; --k) // C - we reverse faster and safer ! we avoid bug if keeps[k] = 0
    {
       size_t keep = keeps[k];
  
        if (keep == 0) {
            S.remove_block(S.left_basis()[k].first,
                           S.right_basis()[k].first);
            U.remove_block(U.left_basis()[k].first,
                           U.right_basis()[k].first);
            V.remove_block(V.left_basis()[k].first,
                           V.right_basis()[k].first);
  // C- idem heev_truncate          --k; // everything gets shifted, to we have to look into the same k again
        } else {
            #ifdef AMBIENT
            ambient::numeric::split(S(S.left_basis()[k].first, S.right_basis()[k].first));
            ambient::numeric::split(U(U.left_basis()[k].first, U.right_basis()[k].first));
            ambient::numeric::split(V(V.left_basis()[k].first, V.right_basis()[k].first));
            #endif

            if (keep >= num_rows(S[k])) continue;
        
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

    delete[] keeps;

    std::size_t bond_dimension = S.left_basis().sum_of_sizes();
    if(verbose){
        maquis::cout << "Sum: " << old_basis.sum_of_sizes() << " -> " << bond_dimension << std::endl;
    }
    
    // MD: for singuler values we care about summing the square of the discraded
    // MD: sum of the discarded values is stored elsewhere
    return truncation_results(bond_dimension, truncated_weight, truncated_fraction, smallest_ev);
}

// TODO: not yet working properly.
template<class Matrix, class DiagMatrix, class SymmGroup>
truncation_results alt_svd_truncate(block_matrix<Matrix, SymmGroup> const & M,
                                    block_matrix<Matrix, SymmGroup> & U,
                                    block_matrix<Matrix, SymmGroup> & V,
                                    block_matrix<DiagMatrix, SymmGroup> & S,
                                    double rel_tol, std::size_t Mmax,
                                    bool verbose = true)
{
    assert( M.left_basis().sum_of_sizes() > 0 && M.right_basis().sum_of_sizes() > 0 );

    block_matrix<Matrix, SymmGroup> t;
    block_matrix<Matrix, SymmGroup> R;
    block_matrix<DiagMatrix, SymmGroup> D;
    
//    maquis::cout << "M:" << std::endl << M;

    
    block_matrix<Matrix, SymmGroup> Mconj = conjugate(M);
    
    gemm(M, transpose(Mconj), t);    
    truncation_results res = heev_truncate(t, U, D, rel_tol, Mmax, verbose);
    
    gemm(transpose(Mconj), U, t);
    qr(t, V, R);
    
    maquis::cout << "S.n_blocks: " << D.n_blocks() << std::endl;
    maquis::cout << "S.trace: " << trace(D) << std::endl;
    maquis::cout << "R.n_blocks: " << R.n_blocks() << std::endl;
    maquis::cout << "R.trace: " << trace(R) << std::endl;
    
//    maquis::cout << "S:" << std::endl << S;
//    maquis::cout << "R:" << std::endl << R;
    
    using std::abs;
    for (int n=0; n<D.n_blocks(); ++n)
        for (int i=0; i<num_rows(D[n]); ++i)
            if ( abs(D[n](i,i) - R[n](i,i)*R[n](i,i)) > 1e-6 )
                maquis::cout << "n=" << n << ", i=" << i << " broken. D=" << D[n](i,i) << ", td=" << R[n](i,i)*R[n](i,i) << std::endl;
    
    S = sqrt(D);
    S /= trace(S);
    V.adjoint_inplace();
    
    maquis::cout << "U:\n" << U.left_basis() << "\n" << U.right_basis() << std::endl;
    maquis::cout << "S:\n" << S.left_basis() << "\n" << S.right_basis() << std::endl;
    maquis::cout << "V:\n" << V.left_basis() << "\n" << V.right_basis() << std::endl;
    
//    assert( td == S );
    return res;
}

template<class Matrix, class DiagMatrix, class SymmGroup>
truncation_results heev_truncate(block_matrix<Matrix, SymmGroup> const & M,
                                 block_matrix<Matrix, SymmGroup> & evecs,
                                 block_matrix<DiagMatrix, SymmGroup> & evals,
                                 double cutoff, std::size_t Mmax,
                                 bool verbose = true)
{
    assert( M.left_basis().sum_of_sizes() > 0 && M.right_basis().sum_of_sizes() > 0 );
    #ifdef AMBIENT
    heev_merged(M, evecs, evals);
    #else
    heev(M, evecs, evals);
    #endif
    Index<SymmGroup> old_basis = evals.left_basis();
    size_t* keeps = new size_t[evals.n_blocks()];
    double truncated_fraction, truncated_weight, smallest_ev;

    estimate_truncation(evals, Mmax, cutoff, keeps, truncated_fraction, truncated_weight, smallest_ev);

    for ( int k = evals.n_blocks() - 1; k >= 0; --k) // C - we reverse faster and safer ! we avoid bug if keeps[k] = 0
    {
        size_t keep = keeps[k];
        
        if (keep == 0) {
            evals.remove_block(evals.left_basis()[k].first,
                               evals.right_basis()[k].first);
            evecs.remove_block(evecs.left_basis()[k].first,
                               evecs.right_basis()[k].first);
//            --k; // everything gets shifted, to we have to look into the same k again
// C - Tim : I reversed the loop because the new version was incompatible with the keeps array, and created a bug when keeps[k]=0.
        } else {
            #ifdef AMBIENT
            ambient::numeric::split(evals(evals.left_basis()[k].first, evals.right_basis()[k].first));
            ambient::numeric::split(evecs(evecs.left_basis()[k].first, evecs.right_basis()[k].first));
            #endif

            if(keep >= num_rows(evals[k])) continue;

            evals.resize_block(evals.left_basis()[k].first,
                               evals.right_basis()[k].first,
                               keep, keep);
            evecs.resize_block(evecs.left_basis()[k].first,
                               evecs.right_basis()[k].first,
                               evecs.left_basis()[k].second,
                               keep);
        }
    }
    delete[] keeps;

    std::size_t bond_dimension = evals.left_basis().sum_of_sizes();
    if(verbose){
        maquis::cout << "Sum: " << old_basis.sum_of_sizes() << " -> " << bond_dimension << std::endl;
    }
    
    // MD: for eigenvalues we care about summing the discraded
    return truncation_results(bond_dimension, truncated_fraction, truncated_fraction, smallest_ev);
}

template<class Matrix, class SymmGroup>
void qr(block_matrix<Matrix, SymmGroup> const& M,
        block_matrix<Matrix, SymmGroup> & Q,
        block_matrix<Matrix, SymmGroup> & R)
{
    /* thin QR in each block */
    Index<SymmGroup> m = M.left_basis(), n = M.right_basis(), k = M.right_basis();
    for (size_t i=0; i<k.size(); ++i)
        k[i].second = std::min(m[i].second,n[i].second);
    
    Q = block_matrix<Matrix, SymmGroup>(m,k);
    R = block_matrix<Matrix, SymmGroup>(k,n);
    std::size_t loop_max = M.n_blocks();
    
    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        qr(M[k], Q[k], R[k]);
    
    assert(Q.right_basis() == R.left_basis());
    assert(Q.reasonable());
    assert(R.reasonable());
}

template<class Matrix, class SymmGroup>
void lq(block_matrix<Matrix, SymmGroup> const& M,
        block_matrix<Matrix, SymmGroup> & L,
        block_matrix<Matrix, SymmGroup> & Q)
{
    /* thin LQ in each block */
    Index<SymmGroup> m = M.left_basis(), n = M.right_basis(), k = M.right_basis();
    for (size_t i=0; i<k.size(); ++i)
        k[i].second = std::min(m[i].second,n[i].second);
    
    L = block_matrix<Matrix, SymmGroup>(m,k);
    Q = block_matrix<Matrix, SymmGroup>(k,n);
    std::size_t loop_max = M.n_blocks();
    
    parallel_for(locale::compact(loop_max), locale k = 0; k < loop_max; ++k)
        lq(M[k], L[k], Q[k]);
    
    assert(Q.left_basis() == L.right_basis());
    assert(Q.reasonable());
    assert(L.reasonable());
}

template<class Matrix, class SymmGroup>
block_matrix<typename maquis::traits::transpose_view<Matrix>::type, SymmGroup> transpose(block_matrix<Matrix, SymmGroup> const & m) 
{ 
    block_matrix<typename maquis::traits::transpose_view<Matrix>::type, SymmGroup> ret; 
    for(size_t k=0; k<m.n_blocks(); ++k) 
        ret.insert_block(transpose(m[k]), m.right_basis()[k].first, m.left_basis()[k].first);
#ifdef AMBIENT_TRACKING
    ambient_track_as(ret, m.label);
#endif
    return ret; 
} 

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> conjugate(block_matrix<Matrix, SymmGroup> m)
{
    m.conjugate_inplace();
    return m;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> adjoint(block_matrix<Matrix, SymmGroup> m)
{
    m.adjoint_inplace();
    return m;
}

template<class Matrix, class SymmGroup>
typename block_matrix<Matrix, SymmGroup>::scalar_type trace(block_matrix<Matrix, SymmGroup> const & m)
{
    return m.trace();
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> adjoin(block_matrix<Matrix, SymmGroup> const & m) // error: it should be adjoin_t_
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
bool is_hermitian(block_matrix<Matrix, SymmGroup> const & m)
{
    bool ret = true;
    for (size_t k=0; ret && k < m.n_blocks(); ++k) {
        if (m.left_basis()[k].second != m.right_basis()[k].second)
            return false;
        else if (m.left_basis()[k].first == m.right_basis()[k].first)
            ret = is_hermitian(m[k]);
        else if (! m.has_block(m.right_basis()[k].first, m.left_basis()[k].first))
            return false;
        else
            ret = ( m[k] == transpose(conj( m(m.right_basis()[k].first, m.left_basis()[k].first) )) );
    }
    return ret;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> sqrt(block_matrix<Matrix, SymmGroup>  m)
{
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        sqrt_inplace(m[k]);

    return m;
}

// Is it really useful?
//template <class Matrix, class SymmGroup, class A>
//block_matrix<Matrix, SymmGroup> exp (block_matrix<Matrix, SymmGroup> const & M, A const & alpha = 1)
//{
//    block_matrix<Matrix, SymmGroup> N, tmp, res;
//    block_matrix<typename alps::numeric::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
//    
//    heev(M, N, S);
//    for (std::size_t k = 0; k < S.n_blocks(); ++k)
//        S[k] = exp(alpha*S[k]);
//    gemm(N, S, tmp);
//    gemm(tmp, adjoint(N), res);
//    
//    return res;
//}

template <class Matrix, class SymmGroup, class A>
block_matrix<Matrix, SymmGroup> op_exp_hermitian(Index<SymmGroup> const & phys,
                                                 block_matrix<Matrix, SymmGroup> M,
                                                 A const & alpha = 1.)
{
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            M(it_c->first, it_c->first) = exp_hermitian(M(it_c->first, it_c->first), alpha);
        else
            M.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return M;
}

namespace detail {
    
    template <class Matrix>
    typename boost::enable_if<boost::is_complex<typename Matrix::value_type>, Matrix>::type
    exp_dispatcher(Matrix const& m, typename Matrix::value_type const& alpha)
    {
        return exp(m, alpha);
    }

    template <class Matrix>
    typename boost::disable_if<boost::is_complex<typename Matrix::value_type>, Matrix>::type
    exp_dispatcher(Matrix const& m, typename Matrix::value_type const& alpha)
    {
        throw std::runtime_error("Exponential of non-hermitian real matrices not implemented!");
        return Matrix();
    }
}

template <class Matrix, class SymmGroup, class A>
block_matrix<Matrix, SymmGroup> op_exp(Index<SymmGroup> const & phys,
                                       block_matrix<Matrix, SymmGroup> M,
                                       A const & alpha = 1.)
{
    for (typename Index<SymmGroup>::const_iterator it_c = phys.begin(); it_c != phys.end(); it_c++)
        if (M.has_block(it_c->first, it_c->first))
            M(it_c->first, it_c->first) = detail::exp_dispatcher(M(it_c->first, it_c->first), alpha);
        else
            M.insert_block(Matrix::identity_matrix(phys.size_of_block(it_c->first)),
                           it_c->first, it_c->first);
    return M;
}

template<class Matrix1, class Matrix2, class SymmGroup>
void op_kron(Index<SymmGroup> const & phys_A,
             Index<SymmGroup> const & phys_B,
             block_matrix<Matrix1, SymmGroup> const & A,
             block_matrix<Matrix1, SymmGroup> const & B,
             block_matrix<Matrix2, SymmGroup> & C)
{
    C = block_matrix<Matrix2, SymmGroup>();

    ProductBasis<SymmGroup> pb_left(phys_A, phys_B);
    ProductBasis<SymmGroup> const& pb_right = pb_left;

    for (int i = 0; i < A.n_blocks(); ++i) {
        for (int j = 0; j < B.n_blocks(); ++j) {
            typename SymmGroup::charge new_right = SymmGroup::fuse(A.right_basis()[i].first, B.right_basis()[j].first);
            typename SymmGroup::charge new_left = SymmGroup::fuse(A.left_basis()[i].first, B.left_basis()[j].first);


            Matrix2 tmp(pb_left.size(A.left_basis()[i].first, B.left_basis()[j].first),
                       pb_right.size(A.right_basis()[i].first, B.right_basis()[j].first),
                       0);

            maquis::dmrg::detail::op_kron(tmp, B[j], A[i],
                                          pb_left(A.left_basis()[i].first, B.left_basis()[j].first),
                                          pb_right(A.right_basis()[i].first, B.right_basis()[j].first),
                                          A.left_basis()[i].second, B.left_basis()[j].second,
                                          A.right_basis()[i].second, B.right_basis()[j].second);

            C.match_and_add_block(tmp, new_left, new_right);
        }
    }
}

template<class Matrix, class SymmGroup>
void op_kron_long(MultiIndex<SymmGroup> const & midx,
                  typename MultiIndex<SymmGroup>::set_id s,
                  block_matrix<Matrix, SymmGroup> const & A,
                  block_matrix<Matrix, SymmGroup> const & B,
                  block_matrix<Matrix, SymmGroup> const & F,
                  std::size_t dist,
                  block_matrix<Matrix, SymmGroup> & C)
{
    assert( midx.size() == 2*(dist+1) );
    C = block_matrix<Matrix, SymmGroup>();
    
    for (size_t run=0; run<2; ++run) {
        
        if (run == 1)
            C.allocate_blocks();
        
        for (index_product_iterator<SymmGroup> it = midx.begin();
             it != midx.end(); ++it)
        {
            bool has_block = A.has_block((*it)[0].first, (*it)[1].first);
            has_block = has_block && B.has_block((*it)[2*dist].first, (*it)[2*dist+1].first);
            for (size_t i=1; has_block && i<dist; ++i)
                has_block = F.has_block((*it)[2*i].first, (*it)[2*i+1].first);
            
            if (!has_block)
                continue;
            
            typename Matrix::value_type val = A((*it)[0], (*it)[1]) * B((*it)[2*dist], (*it)[2*dist+1]);
            for (size_t i=1; i<dist; ++i)
                val *= F((*it)[2*i], (*it)[2*i+1]);
            
            if (val != 0.) {
                typename MultiIndex<SymmGroup>::coord_t coord_l, coord_r;
                boost::tie(coord_l, coord_r) = midx.get_coords(s, *it);
                if (run == 0)
                    C.reserve(coord_l.first, coord_r.first,
                              midx.left_size(s, coord_l.first), midx.right_size(s, coord_r.first));
                else
                    C(coord_l, coord_r) += val;
                
            }
        }
        
    }
    
}

#endif
