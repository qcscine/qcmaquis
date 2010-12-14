#include "mp_tensors/mpstensor.h"

#include "block_matrix/reshapes.h"
#include "block_matrix/block_matrix_algorithms.h"

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>::MPSTensor(Index<SymmGroup> const & sd,
                                        Index<SymmGroup> const & ld,
                                        Index<SymmGroup> const & rd,
                                        bool fillrand)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
//, data_(sd*ld, rd)
, cur_storage(LeftPaired)
, cur_normalization(U)
{
    Index<SymmGroup> lb = sd*ld, rb = rd;
    common_subset(lb, rb);
    
    data_ = block_matrix<Matrix, SymmGroup>(lb, rb);
    
    if (fillrand)
//        data_.generate(drand48);
        data_.generate(rand);
} 

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> MPSTensor<Matrix, SymmGroup>::get_reflected() const
{
    MPSTensor<Matrix, SymmGroup> ret;
    ret.left_i = right_i;
    ret.right_i = left_i;
    ret.phys_i = phys_i;
    
    ret.data_ = transpose(data_);
    
    ret.cur_storage = (cur_storage == LeftPaired ? RightPaired : LeftPaired);
    if (cur_normalization == L)
        ret.cur_normalization = R;
    else if (cur_normalization == R)
        ret.cur_normalization = L;
    else
        ret.cur_normalization = U;
    
    return ret;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & MPSTensor<Matrix, SymmGroup>::site_dim() const
{
    return phys_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & MPSTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & MPSTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isobccompatible(Indicator i) const
{
    return false;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::make_left_paired() const
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
void MPSTensor<Matrix, SymmGroup>::make_right_paired() const
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
    
    
    MPSTensor<Matrix, SymmGroup> rt = get_reflected();
    block_matrix<Matrix, SymmGroup> r = rt.normalize_left(method, multiplied, truncation, bond_dim);
    r = transpose(r);
    *this = rt.get_reflected();
    
    return r;
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_right(block_matrix<Matrix, SymmGroup> const & N)
{
    cur_normalization = U;
    block_matrix<Matrix, SymmGroup> tmp;
    make_left_paired();
    gemm(data_, N, tmp);
    swap(data_, tmp);
    right_i = N.right_basis();
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_left(block_matrix<Matrix, SymmGroup> const & N)
{
    cur_normalization = U;
    block_matrix<Matrix, SymmGroup> tmp;
    make_right_paired();
    gemm(N, data_, tmp);
    swap(data_, tmp);
    left_i = N.left_basis();
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_by_scalar(scalar_type s)
{
    data_ *= s;
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::real_type
MPSTensor<Matrix, SymmGroup>::scalar_norm() const
{
    make_left_paired();
    block_matrix<Matrix, SymmGroup> t;
    gemm(conjugate_transpose(data_), data_, t);
    return sqrt(trace(t));
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::scalar_type
MPSTensor<Matrix, SymmGroup>::scalar_overlap(MPSTensor<Matrix, SymmGroup> const & rhs) const
{
    make_left_paired();
    rhs.make_left_paired();
    block_matrix<Matrix, SymmGroup> t;
    gemm(conjugate_transpose(data_), rhs.data_, t);
    return trace(t);
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

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> &
MPSTensor<Matrix, SymmGroup>::data()
{
    return data_;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::data() const
{
    return data_;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator*=(scalar_type t)
{
    data_ *= t;
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator/=(scalar_type t)
{
    data_ /= t;
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator+=(MPSTensor<Matrix, SymmGroup> const & rhs)
{
    assert(std::equal(left_i.begin(), left_i.end(), rhs.left_i.begin()));
    assert(std::equal(right_i.begin(), right_i.end(), rhs.right_i.begin()));
    assert(std::equal(phys_i.begin(), phys_i.end(), rhs.phys_i.begin()));
    
    make_left_paired();
    rhs.make_left_paired();
    
    data_ += rhs.data_;
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator-=(MPSTensor<Matrix, SymmGroup> const & rhs)
{
    assert(std::equal(left_i.begin(), left_i.end(), rhs.left_i.begin()));
    assert(std::equal(right_i.begin(), right_i.end(), rhs.right_i.begin()));
    assert(std::equal(phys_i.begin(), phys_i.end(), rhs.phys_i.begin()));
    
    make_left_paired();
    rhs.make_left_paired();
    
    data_ -= rhs.data_;
    return *this;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::swap_with(MPSTensor<Matrix, SymmGroup> & b)
{
    swap(this->phys_i, b.phys_i);
    swap(this->left_i, b.left_i);
    swap(this->right_i, b.right_i);
    swap(this->data_, b.data_);
    swap(this->cur_storage, b.cur_storage);
    swap(this->cur_normalization, b.cur_normalization);
}

namespace ietl
{
    template<class Matrix, class SymmGroup, class Generator> void generate(MPSTensor<Matrix, SymmGroup> & m, Generator g)
    {
        m.data().generate(g);
    }
        
    template<class Matrix, class SymmGroup> void generate(MPSTensor<Matrix, SymmGroup> & m, MPSTensor<Matrix, SymmGroup> const & m2)
    {
        m = m2;
    }
    
    template<class Matrix, class SymmGroup> void swap(MPSTensor<Matrix, SymmGroup> & x, MPSTensor<Matrix, SymmGroup> & y)
    {
        x.swap_with(y);
    }
    
    template<class Matrix, class SymmGroup>
    typename MPSTensor<Matrix, SymmGroup>::scalar_type
    dot(MPSTensor<Matrix, SymmGroup> const & x, MPSTensor<Matrix, SymmGroup> const & y)
    {
        return x.scalar_overlap(y);
    }
    
    template<class Matrix, class SymmGroup>
    typename MPSTensor<Matrix, SymmGroup>::scalar_type
    two_norm(MPSTensor<Matrix, SymmGroup> const & x)
    {
        return x.scalar_norm();
    }
}

