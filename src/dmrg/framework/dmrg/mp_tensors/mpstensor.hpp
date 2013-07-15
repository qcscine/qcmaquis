/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/mpstensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/utils/random.hpp"
#include <alps/numeric/isnan.hpp>
#include <alps/numeric/isinf.hpp>

#include "utils/traits.hpp"


// implementation of join functions
#include "dmrg/mp_tensors/joins.hpp"

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>::MPSTensor(Index<SymmGroup> const & sd,
                                        Index<SymmGroup> const & ld,
                                        Index<SymmGroup> const & rd,
                                        bool fillrand,
                                        typename Matrix::value_type val)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
, cur_storage(LeftPaired)
, cur_normalization(Unorm)
{
    Index<SymmGroup> lb = sd*ld, rb = rd;
    common_subset(lb, rb);
    
    // remove blocks from the right index that may not be allowed by the left index
    right_i = rb;
    // remove blocks from the left index that may not be allowed by the right index
    Index<SymmGroup> possible_rp = adjoin(phys_i)*right_i, ltemp = ld;
    common_subset(ltemp, possible_rp);
    left_i = ltemp;
    
    // this is safe, since sd \otimes ld = rd forces the block_matrix to be diagonal
    lb.sort();
    rb.sort();
    left_i.sort();
    right_i.sort();
    
    data_ = block_matrix<Matrix, SymmGroup>(lb, rb);
    
    if (fillrand)
        data_.generate(static_cast<dmrg_random::value_type(*)()>(&dmrg_random::uniform));
    else
        data_.generate(utils::constant<typename Matrix::value_type>(val));
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup>::MPSTensor(Index<SymmGroup> const& sd,
          Index<SymmGroup> const& ld,
          Index<SymmGroup> const& rd,
          block_matrix<Matrix, SymmGroup> const& block,
          MPSStorageLayout layout)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
, data_(block)
, cur_storage(layout)
, cur_normalization(Unorm)
{ }

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::replace_left_paired(block_matrix<Matrix, SymmGroup> const & rhs, Indicator normalization)
{
    make_left_paired();
    
    Index<SymmGroup> new_right_i = rhs.right_basis();
    Index<SymmGroup> possible_left_i = adjoin(phys_i)*new_right_i;
    Index<SymmGroup> old_left_i = left_i;
    
    common_subset(old_left_i, possible_left_i);
    
    right_i = new_right_i;
    left_i = old_left_i;
    
    data_ = rhs;
    
    cur_normalization = normalization;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::replace_right_paired(block_matrix<Matrix, SymmGroup> const & rhs, Indicator normalization)
{
    make_right_paired();
    
    Index<SymmGroup> new_left_i = rhs.left_basis();
    Index<SymmGroup> possible_right_i = phys_i*new_left_i;
    Index<SymmGroup> old_right_i = right_i;
    
    common_subset(old_right_i, possible_right_i);
    
    left_i = new_left_i;
    right_i = old_right_i;
    
    data_ = rhs;
    cur_normalization = normalization;
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
    reshape_right_to_left_new<Matrix>(phys_i, left_i, right_i,
                                      data_, tmp);
    swap(data_, tmp);
    cur_storage = LeftPaired;
    
    assert( weak_equal(right_i, data_.right_basis()) );
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::make_right_paired() const
{   
    if (cur_storage == RightPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    reshape_left_to_right_new<Matrix>(phys_i, left_i, right_i,
                                  data_, tmp);
    swap(data_, tmp);
    cur_storage = RightPaired;
    
    assert( weak_equal(left_i, data_.left_basis()) );
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_left(DecompMethod method,
                                             bool multiplied,
                                             double truncation,
                                             Index<SymmGroup> bond_dim)
{
    if (cur_normalization == Unorm || cur_normalization == Rnorm) {
        if (method == QR) {
            make_left_paired();
            
            block_matrix<Matrix, SymmGroup> Q, R;
            qr(data_, Q, R);
            
            swap(data_, Q);
            right_i = data_.right_basis();
            assert(right_i == R.left_basis());
            
            cur_normalization = Lnorm;
            return R;
        } else {
            make_left_paired();
            
            block_matrix<Matrix, SymmGroup> U, V;
            block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
            
            svd(data_, U, V, S);
            
            right_i = U.right_basis();
            assert(data_.left_basis() == U.left_basis());
            
            swap(data_, U);
            gemm(S, V, U);

            cur_normalization = Lnorm;
            return U;
        }
    }
    return identity_matrix<Matrix>(data_.right_basis());
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup>
MPSTensor<Matrix, SymmGroup>::normalize_right(DecompMethod method,
                                              bool multiplied,
                                              double truncation,
                                              Index<SymmGroup> bond_dim)
{
    if (cur_normalization == Unorm || cur_normalization == Lnorm) {
        if (method == QR) { //enum QR but LQ decomposition
            make_right_paired();
            
            block_matrix<Matrix, SymmGroup> L, Q;
            lq(data_, L, Q);
            
            swap(data_, Q);
            left_i = data_.left_basis();
            assert(left_i == L.right_basis());
            
            cur_normalization = Rnorm;
            return L;
        } else {
            make_right_paired();
            
            block_matrix<Matrix, SymmGroup> U, V;
            block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
            
            svd(data_, U, V, S);
            
            left_i = V.left_basis();
            assert(data_.right_basis() == V.right_basis());
            swap(data_, V);
            
            gemm(U, S, V);
            
            cur_normalization = Rnorm;
            return V;
        }
    }
    return identity_matrix<Matrix>(data_.left_basis());
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::shift_aux_charges(typename SymmGroup::charge diff)
{
    left_i.shift(diff);
    right_i.shift(diff);    
    data_.shift_basis(diff);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_right(block_matrix<Matrix, SymmGroup> const & N)
{
    cur_normalization = Unorm;
    block_matrix<Matrix, SymmGroup> tmp;
    make_left_paired();
    gemm(data_, N, tmp);
    replace_left_paired(tmp);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_from_left(block_matrix<Matrix, SymmGroup> const & N)
{
    cur_normalization = Unorm;
    block_matrix<Matrix, SymmGroup> tmp;
    make_right_paired();
    gemm(N, data_, tmp);
    replace_right_paired(tmp);
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::multiply_by_scalar(const scalar_type& s)
{
    cur_normalization = Unorm;
    *this *= s;
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::divide_by_scalar(const scalar_type& s)
{
    cur_normalization = Unorm;
    *this /= s;
}

template<class Matrix, class SymmGroup>
void
MPSTensor<Matrix, SymmGroup>::conjugate_inplace()
{
    data_ = data_.inplace_conjugate();
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::real_type
MPSTensor<Matrix, SymmGroup>::scalar_norm() const 
{
    return data_.norm();
}

template<class T>
void verbose_assert(T const & a, T const & b)
{
    if (!( a == b) ) {
        maquis::cout << "a: " << a << std::endl;
        maquis::cout << "b: " << b << std::endl;
        assert( a == b );
    }
}

template<class Matrix, class SymmGroup>
typename MPSTensor<Matrix, SymmGroup>::scalar_type
MPSTensor<Matrix, SymmGroup>::scalar_overlap(MPSTensor<Matrix, SymmGroup> const & rhs) const
{
    make_left_paired();
    rhs.make_left_paired();

    // verbose_assert(left_i, rhs.left_i);
    // verbose_assert(right_i, rhs.right_i);
    // verbose_assert(phys_i, rhs.phys_i);
    // verbose_assert(rhs.data_.left_basis(), data_.left_basis());
    // verbose_assert(rhs.data_.right_basis(), data_.right_basis());
    // verbose_assert(rhs.data_.n_blocks(), data_.n_blocks());
    
    // Bela says: this is a workaround for the very rare condition that site_hamil2 removes blocks
    // This shouldn't be necessary, but as of Rev. 1702, is necessary in some cases
    // If I haven't fixed this by the end of Feb 2012, remind me
    Index<SymmGroup> i1 = data_.left_basis(), i2 = rhs.data_.left_basis();
    common_subset(i1, i2);
    std::vector<scalar_type> vt; vt.reserve(i1.size());

    for (std::size_t b = 0; b < i1.size(); ++b) {
        typename SymmGroup::charge c = i1[b].first;
        assert( data_.has_block(c,c) && rhs.data_.has_block(c,c) );
        vt.push_back(overlap(data_(c,c), rhs.data_(c,c)));
    } // should be reformulated in terms of reduction (todo: Matthias, 30.04.12 / scalar-value types)
    return maquis::accumulate(vt.begin(), vt.end(), scalar_type(0.));
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
bool MPSTensor<Matrix, SymmGroup>::isleftnormalized(bool test) const
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == Lnorm;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isrightnormalized(bool test) const
{
    if (test)
        throw std::runtime_error("Not implemented!");
    else
        return cur_normalization == Rnorm;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::isnormalized(bool test) const
{
    if (isleftnormalized(test) || isrightnormalized(test))
        return true;
    else
        return false;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> &
MPSTensor<Matrix, SymmGroup>::data()
{
    cur_normalization = Unorm;
    return data_;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::data() const
{
    return const_data();
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::const_data() const
{
    return data_;
}


template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator*=(const scalar_type& t)
{
    cur_normalization = Unorm;
    data_ *= t;
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator/=(const scalar_type& t)
{
    cur_normalization = Unorm;
    data_ /= t;
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator+=(MPSTensor<Matrix, SymmGroup> const & rhs)
{
    assert( weak_equal(left_i, rhs.left_i) );
    assert( weak_equal(right_i, rhs.right_i) );
    assert( phys_i == rhs.phys_i );
    
    make_left_paired();
    rhs.make_left_paired();
    
    cur_normalization = Unorm;

    for (std::size_t i = 0; i < data_.n_blocks(); ++i)
    {
        typename SymmGroup::charge lc = data_.left_basis()[i].first, rc = data_.right_basis()[i].first;
        std::size_t matched_block = rhs.data_.find_block(lc,rc);
        if (matched_block < rhs.data_.n_blocks()) {
            data_[i] += rhs.data_[matched_block];
        }
    }
    
    return *this;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> const &
MPSTensor<Matrix, SymmGroup>::operator-=(MPSTensor<Matrix, SymmGroup> const & rhs)
{
    assert( weak_equal(left_i, rhs.left_i) );
    assert( weak_equal(right_i, rhs.right_i) );
    assert( phys_i == rhs.phys_i );
    
    make_left_paired();
    rhs.make_left_paired();
    
    cur_normalization = Unorm;
    
    for (std::size_t i = 0; i < data_.n_blocks(); ++i)
    {
        typename SymmGroup::charge lc = data_.left_basis()[i].first, rc = data_.right_basis()[i].first;
        if (rhs.data_.has_block(lc,rc)) {
            data_[i] -= rhs.data_(lc,rc);
        }
    }
    
    return *this;
}

template<class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::swap_with(MPSTensor<Matrix, SymmGroup> & b)
{
    using std::swap;
    swap(this->phys_i, b.phys_i);
    swap(this->left_i, b.left_i);
    swap(this->right_i, b.right_i);
    swap(this->data_, b.data_);
    swap(this->cur_storage, b.cur_storage);
    swap(this->cur_normalization, b.cur_normalization);
}

template<class Matrix, class SymmGroup>
template<class Archive>
void MPSTensor<Matrix, SymmGroup>::load(Archive & ar)
{
    data_.clear();
    make_left_paired();
    ar["phys_i"] >> phys_i;
    ar["left_i"] >> left_i;
    ar["right_i"] >> right_i;
    ar["data_"] >> data_;
    cur_normalization = Unorm;
}

template<class Matrix, class SymmGroup>
template<class Archive>
void MPSTensor<Matrix, SymmGroup>::save(Archive & ar) const
{
    make_left_paired();
    ar["phys_i"] << phys_i;
    ar["left_i"] << left_i;
    ar["right_i"] << right_i;
    ar["data_"] << data_;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::reasonable() const
{
    {
        make_left_paired();
        if ( !weak_equal(right_i, data_.right_basis()) )
            throw std::runtime_error("right basis is wrong");
        
//        maquis::cout << "** reasonable left_paired **" << std::endl;
//        maquis::cout << "reasonable::left_i: " << left_i << std::endl;
//        maquis::cout << "reasonable::right_i: " << right_i << std::endl;
//        maquis::cout << "reasonable::data_:" << std::endl << data_ << std::endl;
        make_right_paired();
        if ( !weak_equal(left_i, data_.left_basis()) )
            throw std::runtime_error("left basis is wrong");
        
//        maquis::cout << "** reasonable right_paired **" << std::endl;
//        maquis::cout << "reasonable::left_i: " << left_i << std::endl;
//        maquis::cout << "reasonable::right_i: " << right_i << std::endl;
//        maquis::cout << "reasonable::data_:" << std::endl << data_ << std::endl;
    }
    
    {
        for (std::size_t i = 0; i < data_.n_blocks(); ++i)
        {
            if (data_.left_basis()[i].first != data_.right_basis()[i].first)
                throw std::runtime_error("particle number is wrong");
        }
    }
    return true;
}

template<class Matrix, class SymmGroup>
bool MPSTensor<Matrix, SymmGroup>::num_check() const
{
        for (std::size_t k = 0; k < data_.n_blocks(); ++k)
        {
            for (size_t i = 0; i<num_rows(data_[k]); ++i)
                for (size_t j = 0; j<num_cols(data_[k]); ++j)
                {
                    if ( alps::numeric::isnan(data_[k](i,j)) )
                        throw std::runtime_error("NaN found!");
                    if ( alps::numeric::isinf(data_[k](i,j)) )
                        throw std::runtime_error("INF found!");
                }
        }
    return true;
}

template <class Matrix, class SymmGroup>
void MPSTensor<Matrix, SymmGroup>::check_equal (MPSTensor<Matrix, SymmGroup> const & rhs) const
{
    std::string error;
    // Indexes
    
    try {
        reasonable();
    } catch (std::exception & e) {
        error += "I'm unreasonable. " + std::string(e.what());
    }
    try {
        rhs.reasonable();
    } catch (std::exception & e) {
        error += "He's unreasonable. " + std::string(e.what());
    }
    
    // Storage
//    if (cur_storage != rhs.cur_storage)
//        error += "Storage doesn't match. ";

    make_left_paired();
    rhs.make_left_paired();
    
    // Data
    if (data_.n_blocks() != rhs.data_.n_blocks())
        error += "n_blocks doesn't match. ";
    else {
        for (int b=0; b<data_.n_blocks(); ++b) {
            if (data_[b].num_cols() != rhs.data_[b].num_cols() || data_[b].num_rows() != rhs.data_[b].num_rows())
                error += "Size of block doesn't match. ";
            for (int i=0; i<data_[b].num_rows() && error.empty(); ++i)
                for (int j=0; j<data_[b].num_cols() && error.empty(); ++j)
                    if (data_[b](i,j) != rhs.data_[b](i,j))
                        error += "Data doesn't match. ";
        }
    }
    // Finalize
    if (!error.empty())
        throw std::runtime_error(error);
        
}

template<class Matrix, class SymmGroup>
std::size_t MPSTensor<Matrix, SymmGroup>::num_elements() const
{
    return data_.num_elements();
}

