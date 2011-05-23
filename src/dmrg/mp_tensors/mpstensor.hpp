/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "mp_tensors/mpstensor.h"

#include "mp_tensors/reshapes.h"
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
, cur_normalization(Unorm)
{
    Index<SymmGroup> lb = sd*ld, rb = rd;
    common_subset(lb, rb);
    
    data_ = block_matrix<Matrix, SymmGroup>(lb, rb);
    
    if (fillrand) {
        data_.generate(drand48);
//        for (int i = 0; i < data_.n_blocks(); ++i)
//            data_[i](0,0) = 1;
//            for (int k = 0; k < std::min(num_rows(data_[i]), num_cols(data_[i])); ++k)
//                data_[i](k,k) = 1;
    } else
        data_.generate(utils::constant<typename Matrix::value_type>(0));
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
        throw std::runtime_error("Not implemented!");
        make_left_paired();
        
        block_matrix<Matrix, SymmGroup> q, r;
        qr(data_, q, r);
        swap(data_, q);
        
        cur_normalization = Lnorm;
        return r;
    } else {
        make_left_paired();
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;
        
        svd(data_, U, V, S);
        
        right_i = U.right_basis();
        assert(data_.left_basis() == U.left_basis());
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
    if (method == QR) {
        throw std::runtime_error("Not implemented.");
        
    } else {
        make_right_paired();
        
        block_matrix<Matrix, SymmGroup> U, V;
        block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S;

        svd(data_, U, V, S);
        
        left_i = V.left_basis();
        assert(data_.right_basis() == V.right_basis());
        swap(data_, V);
        
        gemm(U, S, V);
        return V;
    }
    
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_right(block_matrix<Matrix, SymmGroup> const & N)
{
    cur_normalization = Unorm;
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
    cur_normalization = Unorm;
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
    static Timer timer("scalar_norm");
    timer.begin();

//    make_left_paired();
//    block_matrix<Matrix, SymmGroup> t;
//    pgemm(conjugate_transpose(data_), data_, t);
//    real_type r = sqrt(trace(t));

    using utils::conj;
//    printf("<<< Scalar norm playouts\n");
    scalar_type ret = 0;
    for (std::size_t b = 0; b < data_.n_blocks(); ++b)
        for (std::size_t c = 0; c < num_cols(data_[b]); ++c)
            for (std::size_t r = 0; r < num_rows(data_[b]); ++r)
                ret += conj(data_[b](r,c)) * data_[b](r,c);
    
//    printf("Scalar norm playouts >>>\n");
    timer.end();
    return sqrt(ret);
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::scalar_type
MPSTensor<Matrix, SymmGroup>::scalar_overlap(MPSTensor<Matrix, SymmGroup> const & rhs) const
{
    static Timer timer("scalar_overlap");
    timer.begin();

    make_left_paired();
    rhs.make_left_paired();
    
    using utils::conj;
    
    scalar_type ret = 0;
    for (std::size_t b = 0; b < data_.n_blocks(); ++b)
        for (std::size_t c = 0; c < num_cols(data_[b]); ++c)
            for (std::size_t r = 0; r < num_rows(data_[b]); ++r)
                ret += conj(data_[b](r,c)) * rhs.data_[b](r,c);

    timer.end();
    return ret;
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, MPSTensor<Matrix, SymmGroup> const & mps)
{
    os << "Physical space: " << mps.phys_i << std::endl;
    os << "Left space: " << mps.left_i << std::endl;
    os << "Right space: " << mps.right_i << std::endl;
    os << mps.data_;
    return os;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isleftnormalized(bool test)
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == Lnorm;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isrightnormalized(bool test)
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == Rnorm;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isnormalized(bool test)
{
    if (isleftnormalized(test) || isrightnormalized(test))
        return true;
    else {
        cur_normalization = Unorm;
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

#ifdef HAVE_ALPS_HDF5

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::serialize(alps::hdf5::iarchive & ar)
{
    make_left_paired();
    ar >> alps::make_pvp("phys_i", phys_i);
    ar >> alps::make_pvp("left_i", left_i);
    ar >> alps::make_pvp("right_i", right_i);
    ar >> alps::make_pvp("data_", data_);
    cur_normalization = Unorm;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::serialize(alps::hdf5::oarchive & ar) const
{
    make_left_paired();
    ar << alps::make_pvp("phys_i", phys_i);
    ar << alps::make_pvp("left_i", left_i);
    ar << alps::make_pvp("right_i", right_i);
    ar << alps::make_pvp("data_", data_);
}

#endif

