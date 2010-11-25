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
{ } 

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
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_left(DecompMethod method,
                                             bool multiplied,
                                             double truncation,
                                             Index<SymmGroup> bond_dim)
{
    if (method != QR)
        throw std::runtime_error("Not implemented!");
    
    make_left_paired();
    
    block_matrix<Matrix, SymmGroup> q, r;
    qr(data_, q, r);
    swap(data_, q);
    
    cur_normalization = L;
    return r;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_right(DecompMethod method,
                                              bool multiplied,
                                              double truncation,
                                              Index<SymmGroup> bond_dim)
{
    if (method != QR)
        throw std::runtime_error("Not implemented!");
    
    make_right_paired();
    
    // optimization idea: also implement rq decomposition
    
    transpose(data_);
    
    block_matrix<Matrix, SymmGroup> q, r;
    qr(data_, q, r);
    swap(data_, q);
    
    transpose(data_);
    transpose(r);
    
    cur_normalization = R;
    return r;
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_right(block_matrix<Matrix, SymmGroup> & N)
{
    block_matrix<Matrix, SymmGroup> tmp;
    make_left_paired();
    gemm(data_, N, tmp);
    swap(data_, tmp);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_left(block_matrix<Matrix, SymmGroup> & N)
{
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
    block_matrix<Matrix, SymmGroup> cp = adjoin(data_), t;
    gemm(cp, data_, t);
    return trace(t);
}
