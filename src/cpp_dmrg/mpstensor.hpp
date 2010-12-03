#include "mpstensor.h"

#include "reshapes.h"
#include "block_matrix_algorithms.h"

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>::MPSTensor(Index<SymmGroup> const & sd,
                                        Index<SymmGroup> const & ld,
                                        Index<SymmGroup> const & rd)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
, data_(sd*ld, rd)
, cur_storage(LeftPaired)
, cur_normalization(U)
{
    data_.fill(drand48);
} 

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::reflect()
{
    swap(left_i, right_i);
    data_ = transpose(data_);
    cur_storage = (cur_storage == LeftPaired ? RightPaired : LeftPaired);
    if (cur_normalization == L)
        cur_normalization = R;
    else if (cur_normalization == R)
        cur_normalization = L;
    // if U, nothing changes
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> MPSTensor<Matrix, SymmGroup>::site_dim() const
{
    return phys_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> MPSTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> MPSTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isobccompatible(Indicator i) const
{
    return false;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::make_left_paired()
{
    if (cur_storage == LeftPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    reshape_right_to_left<Matrix>(phys_i, left_i, right_i,
                                  data_, tmp);
    swap(data_, tmp);
    cur_storage = LeftPaired;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::make_right_paired()
{
    if (cur_storage == RightPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    reshape_left_to_right<Matrix>(phys_i, left_i, right_i,
                                  data_, tmp);
    swap(data_, tmp);
    cur_storage = RightPaired;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_left(DecompMethod method,
                                             bool multiplied,
                                             double truncation,
                                             Index<SymmGroup> bond_dim)
{
    if (method == QR) {
        make_left_paired();
        
        block_matrix<Matrix, SymmGroup> q, r;
        qr(data_, q, r);
        swap(data_, q);
        
        cur_normalization = L;
        return r;
    } else {
        make_left_paired();
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
        
        svd(data_, U, V, S);
        
        swap(data_, U);
        gemm(S, V, U);
        return U;
    }
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_right(DecompMethod method,
                                              bool multiplied,
                                              double truncation,
                                              Index<SymmGroup> bond_dim)
{
    // This is more expensive than a full implementation
    // due to the additional transposes - but I want to reduce code for now.
    reflect();
    block_matrix<Matrix, SymmGroup> r = normalize_left(method, multiplied, truncation, bond_dim);
    r = transpose(r);
    reflect();
    
    return r;
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_right(block_matrix<Matrix, SymmGroup> & N)
{
    cur_normalization = U;
    block_matrix<Matrix, SymmGroup> tmp;
    make_left_paired();
    gemm(data_, N, tmp);
    swap(data_, tmp);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_left(block_matrix<Matrix, SymmGroup> & N)
{
    cur_normalization = U;
    block_matrix<Matrix, SymmGroup> tmp;
    make_right_paired();
    gemm(N, data_, tmp);
    swap(data_, tmp);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_by_scalar(scalar_type s)
{
    data_ *= s;
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::real_type
MPSTensor<Matrix, SymmGroup>::scalar_norm()
{
    make_left_paired();
    block_matrix<Matrix, SymmGroup> cp = conjugate_transpose(data_), t;
    gemm(cp, data_, t);
    return sqrt(trace(t));
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, MPSTensor<Matrix, SymmGroup> const & mps)
{
    os << "Physical space: " << mps.phys_i << endl;
    os << "Left space: " << mps.left_i << endl;
    os << "Right space: " << mps.right_i << endl;
    os << mps.data_;
    return os;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isleftnormalized(bool test)
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == L;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isrightnormalized(bool test)
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == R;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isnormalized(bool test)
{
    if (isleftnormalized(test) || isrightnormalized(test))
        return true;
    else {
        cur_normalization = U;
        return false;
    }
}
