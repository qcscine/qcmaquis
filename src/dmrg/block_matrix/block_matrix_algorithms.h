#ifndef BLOCK_MATRIX_ALGORITHMS_H
#define BLOCK_MATRIX_ALGORITHMS_H

#include "block_matrix/block_matrix.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"

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
        
        C.insert_block(boost::tuples::make_tuple(Matrix(num_rows(A[k]), num_columns(B[matched_block])),
                                                 A.left_basis()[k].first, B.right_basis()[matched_block].first));
        gemm(A[k], B[matched_block], C[C.left_basis().position(A.left_basis()[k].first)]);
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
    
    for (std::size_t k = 0; k < M.n_blocks(); ++k)
        svd(M[k], U[k], V[k], S[k]);
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd(block_matrix<Matrix, SymmGroup> const & M,
         block_matrix<Matrix, SymmGroup> & U,
         block_matrix<Matrix, SymmGroup> & V,
         block_matrix<DiagMatrix, SymmGroup> & S,
         double rel_tol)
{
    svd(M, U, V, S);
    
    /* to be done:
     Given the full SVD in each block (above), remove all singular values and corresponding rows/cols
     where the singular value is < rel_tol*max(S), where the maximum is taken over all blocks.
     Be careful to update the Index descriptions in the matrices to reflect the reduced block sizes
     (remove_rows/remove_columns methods for that)
     */
    
    Index<SymmGroup> old_basis = S.left_basis();
    
    double Smax = 0;
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
        Smax = std::max(Smax, *std::max_element(S[k].elements().first, S[k].elements().second));
    double Scut = rel_tol * Smax;
    
    for (std::size_t k = 0; k < S.n_blocks(); ++k)
    {
        std::size_t keep = std::find_if(S[k].elements().first, S[k].elements().second,
                                        boost::lambda::_1 < Scut)-S[k].elements().first;
        if (keep >= num_rows(S[k]))
            continue;
        if (keep == 0)
            keep = 1;
        
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
    
    if (! (old_basis == S.left_basis()) ) {
        cout << "SVD performed a truncation: " << endl;
//        cout << old_basis << endl << S.left_basis() << endl;
        cout << old_basis.sum_of_sizes() << " -> " << S.left_basis().sum_of_sizes() << endl;
    }
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
        ret.insert_block(boost::tuples::make_tuple(m[k],
                                                   -m.left_basis()[k].first,
                                                   -m.right_basis()[k].first));
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
        ret[k] = blas::identity_matrix<Matrix>(size[k].second);
    return ret;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> sqrt(block_matrix<Matrix, SymmGroup> m)
{
    for (std::size_t k = 0; k < m.n_blocks(); ++k)
        m[k] = sqrt(m[k]);
    return m;
}

#endif
