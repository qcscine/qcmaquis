#ifndef BLOCK_MATRIX_ALGORITHMS_H
#define BLOCK_MATRIX_ALGORITHMS_H

#include "block_matrix.h"
#include "matrix_algorithms.hpp"

// some example functions
template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
void gemm(block_matrix<Matrix1, SymmGroup> & A,
          block_matrix<Matrix2, SymmGroup> & B,
          block_matrix<Matrix3, SymmGroup> & C)
{
    C = block_matrix<Matrix3, SymmGroup>(A.left_basis(), B.right_basis());
    
    // Some checks
    // We should discuss whether we want to use asserts or something that doesn't disappear upon NDEBUG
    assert(A.n_blocks() == B.n_blocks() && B.n_blocks() == C.n_blocks());
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        // is the charge of the column of the first the same as the row of the second?
        assert(A.right_basis()[k].first == B.left_basis()[k].first);
        // is the dimension of the sectors the same?
        assert(A.right_basis()[k].second == B.left_basis()[k].second);
    }
    
    for (std::size_t k = 0; k < A.n_blocks(); ++k)
        //        C[k] = A[k] * B[k];
        gemm(A[k], B[k], C[k]);
}

template<class Matrix, class DiagMatrix, class SymmGroup>
void svd(block_matrix<Matrix, SymmGroup> & M,
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
typename Matrix::value_type trace(block_matrix<Matrix, SymmGroup> m)
{
    return m.trace();
}

template<class Matrix, class SymmGroup, class Generator>
void fill_with_random(block_matrix<Matrix, SymmGroup> & m, Generator & g)
{
    m.fill_with_random(g);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> identity_matrix(Index<SymmGroup> const & size)
{
    block_matrix<Matrix, SymmGroup> ret(size, size);
    for (std::size_t k = 0; k < ret.n_blocks(); ++k)
        ret[k] = blas::identity_matrix<Matrix>(size[k].second);
    return ret;
}

#endif
