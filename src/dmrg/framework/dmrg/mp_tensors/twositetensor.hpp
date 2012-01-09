/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpstensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/generic_reshape.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/utils/random.hpp"
#include <alps/numeric/real.hpp>

#include <vector>
#include <utility>

#include "ts_reshape.h"

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup>::TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & mps1,
                                                MPSTensor<Matrix, SymmGroup> const & mps2)
: phys_i( mps1.site_dim()*mps1.site_dim() )
, phys_i_orig( mps1.site_dim() )
, left_i( mps1.row_dim() )
, right_i( mps2.col_dim() )
, cur_storage(TSBothPaired)
, cur_normalization(Unorm)
{
    assert(mps1.site_dim() == mps2.site_dim());
    
    mps1.make_left_paired();
    mps2.make_right_paired();
    gemm(mps1.data(), mps2.data(), data_);
 
}


template<class Matrix, class SymmGroup>
Index<SymmGroup> const & TwoSiteTensor<Matrix, SymmGroup>::site_dim() const
{
    return phys_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & TwoSiteTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
Index<SymmGroup> const & TwoSiteTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::make_left_paired() const
{
    if (cur_storage == TSLeftPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    if (cur_storage == TSBothPaired) {
        ts_reshape::reshape_both_to_left<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    } else {
	// direct left to right reshape should not be needed
	make_both_paired();
        ts_reshape::reshape_both_to_left<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    }
    
    swap(data_, tmp);
    cur_storage = TSLeftPaired;
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::make_both_paired() const
{
    if (cur_storage == TSBothPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    if (cur_storage == TSRightPaired) {
        ts_reshape::reshape_right_to_both<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    }
    else {
        ts_reshape::reshape_left_to_both<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    }
    
    swap(data_, tmp);
    cur_storage = TSBothPaired;
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::make_right_paired() const
{   
    if (cur_storage == TSRightPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    if (cur_storage == TSBothPaired)
        ts_reshape::reshape_both_to_right<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    else {
	// direct left to right reshape should not be needed
	make_both_paired();
        ts_reshape::reshape_both_to_right<Matrix, SymmGroup>(phys_i_orig, left_i, right_i, data_, tmp);
    }
    
    swap(data_, tmp);
    cur_storage = TSRightPaired;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> TwoSiteTensor<Matrix, SymmGroup>::make_mps() const
{
    MPSTensor<Matrix, SymmGroup> tmp(phys_i, left_i, right_i, false, 0);
    make_left_paired();
    tmp.data() = data_;
    return tmp;
}

template<class Matrix, class SymmGroup>
std::pair<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> >
TwoSiteTensor<Matrix, SymmGroup>::split_mps_l2r(std::size_t Mmax, double cutoff, Logger * iter_log) const
{
    make_both_paired();
    
    typedef typename maquis::types::associated_real_diagonal_matrix<Matrix>::type dmt;
    block_matrix<Matrix, SymmGroup> u, v;
    block_matrix<dmt, SymmGroup> s;
    
    svd_truncate(data_, u, v, s, cutoff, Mmax, true, iter_log);
    
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_orig, left_i, right_i, false, 0),
                                 mps_tensor2(phys_i_orig, left_i, right_i, false, 0);
    mps_tensor1.make_left_paired();
    mps_tensor2.make_right_paired();

    mps_tensor1.data() = u;
    mps_tensor1.right_i = u.right_basis();

    gemm(s, v, u);
    mps_tensor2.data() = u;
    mps_tensor2.left_i = u.left_basis();
    
    return std::make_pair(mps_tensor1, mps_tensor2);
}

template<class Matrix, class SymmGroup>
std::pair<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> >
TwoSiteTensor<Matrix, SymmGroup>::split_mps_r2l(std::size_t Mmax, double cutoff, Logger * iter_log) const
{
    make_both_paired();
    
    typedef typename maquis::types::associated_real_diagonal_matrix<Matrix>::type dmt;
    block_matrix<Matrix, SymmGroup> u, v;
    block_matrix<dmt, SymmGroup> s;
    
    svd_truncate(data_, u, v, s, cutoff, Mmax, true, iter_log);
    
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_orig, left_i, right_i, false, 0),
                                 mps_tensor2(phys_i_orig, left_i, right_i, false, 0);
    mps_tensor1.make_left_paired();
    mps_tensor2.make_right_paired();

    mps_tensor2.data() = v;
    mps_tensor2.left_i = v.left_basis();

    gemm(u, s, v);
    mps_tensor1.data() = v;
    mps_tensor1.right_i = v.right_basis();
    
    return std::make_pair(mps_tensor1, mps_tensor2);
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, TwoSiteTensor<Matrix, SymmGroup> const & mps)
{
    os << "Physical space: " << mps.phys_i << std::endl;
    os << "Left space: " << mps.left_i << std::endl;
    os << "Right space: " << mps.right_i << std::endl;
    os << mps.data_;
    return os;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> &
TwoSiteTensor<Matrix, SymmGroup>::data()
{
    return data_;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const &
TwoSiteTensor<Matrix, SymmGroup>::data() const
{
    return data_;
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::swap_with(TwoSiteTensor<Matrix, SymmGroup> & b)
{
    swap(this->phys_i, b.phys_i);
    swap(this->phys_i_orig, b.phys_i_orig);
    swap(this->left_i, b.left_i);
    swap(this->right_i, b.right_i);
    swap(this->data_, b.data_);
    swap(this->cur_storage, b.cur_storage);
    swap(this->cur_normalization, b.cur_normalization);
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> & TwoSiteTensor<Matrix, SymmGroup>::operator << (MPSTensor<Matrix, SymmGroup> const & rhs)
{
    this->make_left_paired();
    rhs.make_left_paired();

    // Precondition: rhs.data() and this->data() have same shape if both are left_paired
    assert( rhs.row_dim() == this->row_dim() &&
		 rhs.col_dim() == this->col_dim() &&
		 rhs.site_dim() == this->site_dim() );
    assert( rhs.data().left_basis() == this->data().left_basis() &&
		 rhs.data().right_basis() == this->data().right_basis() );
    
    this->data() = rhs.data();

    return *this;

}


