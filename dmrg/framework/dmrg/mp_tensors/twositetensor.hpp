/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/mpstensor.h"

#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/utils/random.hpp"
#include <alps/numeric/real.hpp>

#include <vector>
#include <utility>

#include "ts_reshape.h"

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup>::TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & mps1,
                                                MPSTensor<Matrix, SymmGroup> const & mps2)
: phys_i( mps1.site_dim()*mps2.site_dim() )
, phys_i_left( mps1.site_dim() )
, phys_i_right( mps2.site_dim() )
, left_i( mps1.row_dim() )
, right_i( mps2.col_dim() )
, cur_storage(TSBothPaired)
, cur_normalization(Unorm)
{
    mps1.make_left_paired();
    mps2.make_right_paired();
    gemm(mps1.data(), mps2.data(), data_, parallel::scheduler_balanced(mps1.data()));
 
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
Index<SymmGroup> const & TwoSiteTensor<Matrix, SymmGroup>::local_site_dim(short d) const
{
    return (d == 1) ? phys_i_right : phys_i_left;
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::make_left_paired() const
{
    if (cur_storage == TSLeftPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    if (cur_storage == TSBothPaired) {
        ts_reshape::reshape_both_to_left<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    } else {
	    // direct left to right reshape should not be needed
	    make_both_paired();
        ts_reshape::reshape_both_to_left<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    }
    
    swap(data_, tmp);
    cur_storage = TSLeftPaired;
    
    // assert( right_i == data_.right_basis() );
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::make_both_paired() const
{
    if (cur_storage == TSBothPaired)
        return;
    
    block_matrix<Matrix, SymmGroup> tmp;
    if (cur_storage == TSRightPaired) {
        ts_reshape::reshape_right_to_both<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    }
    else {
        ts_reshape::reshape_left_to_both<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
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
        ts_reshape::reshape_both_to_right<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    else {
        // direct left to right reshape should not be needed
	    make_both_paired();
        ts_reshape::reshape_both_to_right<Matrix, SymmGroup>(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    }
    
    swap(data_, tmp);
    cur_storage = TSRightPaired;
    
    // assert( left_i == data_.left_basis() );
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> TwoSiteTensor<Matrix, SymmGroup>::make_mps() const
{
    make_right_paired();
    return MPSTensor<Matrix, SymmGroup>(phys_i, left_i, right_i, data_, RightPaired);
}

template<class Matrix, class SymmGroup>
boost::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::split_mps_l2r(std::size_t Mmax, double cutoff) const
{
    make_both_paired();
    
    typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
    block_matrix<Matrix, SymmGroup> u, v;
    block_matrix<dmt, SymmGroup> s;
    
    truncation_results trunc = svd_truncate(data_, u, v, s, cutoff, Mmax, true);
    
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), u, LeftPaired);
    assert( mps_tensor1.reasonable() );
    gemm(s, v, u);
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, u.left_basis(), right_i, u, RightPaired);
    assert( mps_tensor2.reasonable() );
    
    return boost::make_tuple(mps_tensor1, mps_tensor2, trunc);
}

template<class Matrix, class SymmGroup>
boost::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::split_mps_r2l(std::size_t Mmax, double cutoff) const
{
    make_both_paired();
    
    typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
    block_matrix<Matrix, SymmGroup> u, v;
    block_matrix<dmt, SymmGroup> s;
    
    truncation_results trunc = svd_truncate(data_, u, v, s, cutoff, Mmax, true);
    
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, v.left_basis(), right_i, v, RightPaired);
    
    gemm(u, s, v);
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), v, LeftPaired);
    
    return boost::make_tuple(mps_tensor1, mps_tensor2, trunc);
}

template<class Matrix, class SymmGroup>
boost::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::predict_split_l2r(std::size_t Mmax, double cutoff, double alpha, Boundary<Matrix, SymmGroup> const& left, MPOTensor<Matrix, SymmGroup> const& mpo)
{
    make_both_paired();
    
    /// build reduced density matrix (with left index open)
    block_matrix<Matrix, SymmGroup> dm;
    gemm(data_, transpose(conjugate(data_)), dm, parallel::scheduler_balanced(data_));
    
    /// state prediction
    if (alpha != 0.) {
        maquis::cout << "Growing, alpha = " << alpha << std::endl;
        Index<SymmGroup> right_phys_i = adjoin(phys_i_right) * right_i;
        MPSTensor<Matrix, SymmGroup> tmp(phys_i_left, left_i, right_phys_i, data_, LeftPaired);
        
        Boundary<Matrix, SymmGroup> half_dm = contraction::Engine<Matrix, Matrix, SymmGroup>::left_boundary_tensor_mpo(tmp, left, mpo);
        tmp = MPSTensor<Matrix, SymmGroup>();
        
        omp_for(std::size_t b, parallel::range<std::size_t>(0,half_dm.aux_dim()), {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(half_dm[b], transpose(conjugate(half_dm[b])), tdm);
            tdm *= alpha;
            swap(tdm, half_dm[b]);
        });
        
        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
            block_matrix<Matrix, SymmGroup> const& tdm = half_dm[b];
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (data_.basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
                    dm.reserve(tdm.basis().left_charge(k), tdm.basis().right_charge(k),
                               num_rows(tdm[k]), num_cols(tdm[k]));
            }
        }
        dm.allocate_blocks();
        
        omp_for(std::size_t k, parallel::range<std::size_t>(0,dm.n_blocks()), {
            for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
                std::size_t match = half_dm[b].find_block(dm.basis().left_charge(k), dm.basis().right_charge(k));
                if (match < half_dm[b].n_blocks())
                    dm[k] += (half_dm[b][match]);
            }
        });
    }
    assert( weak_equal(dm.left_basis(), data_.left_basis()) );
    
    /// truncation
    block_matrix<Matrix, SymmGroup> U;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
    truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax);
    dm = block_matrix<Matrix, SymmGroup>();

    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, U.right_basis(), U, LeftPaired);
    assert( mps_tensor1.reasonable() );
    
//    block_matrix<Matrix, SymmGroup> t;
//    gemm(U, S, t);
//    maquis::cout << "t = " << t << std::endl;

    block_matrix<Matrix, SymmGroup> V;
    gemm(transpose(conjugate(U)), data_, V);
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, V.left_basis(), right_i, V, RightPaired);
    assert( mps_tensor2.reasonable() );
    
    return boost::make_tuple(mps_tensor1, mps_tensor2, trunc);
}


template<class Matrix, class SymmGroup>
boost::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::predict_split_r2l(std::size_t Mmax, double cutoff, double alpha, Boundary<Matrix, SymmGroup> const& right, MPOTensor<Matrix, SymmGroup> const& mpo)
{
    make_both_paired();
    
    /// build reduced density matrix (with right index open)
    block_matrix<Matrix, SymmGroup> dm;
    gemm(transpose(conjugate(data_)), data_, dm, parallel::scheduler_balanced(data_));
    
    /// state prediction
    if (alpha != 0.) {
        maquis::cout << "Growing, alpha = " << alpha << std::endl;
        Index<SymmGroup> left_phys_i = phys_i_left * left_i;
        MPSTensor<Matrix, SymmGroup> tmp(phys_i_right, left_phys_i, right_i, data_, RightPaired);
        
        Boundary<Matrix, SymmGroup> half_dm = contraction::Engine<Matrix, Matrix, SymmGroup>::right_boundary_tensor_mpo(tmp, right, mpo);
        tmp = MPSTensor<Matrix, SymmGroup>();
        
        omp_for(std::size_t b, parallel::range<std::size_t>(0,half_dm.aux_dim()), {
            block_matrix<Matrix, SymmGroup> tdm;
            gemm(transpose(conjugate(half_dm[b])), half_dm[b], tdm);
            tdm *= alpha;
            swap(tdm, half_dm[b]);
        });

        for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
            block_matrix<Matrix, SymmGroup> const& tdm = half_dm[b];
            for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
                if (data_.basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
                    dm.reserve(tdm.basis().left_charge(k), tdm.basis().right_charge(k),
                               num_rows(tdm[k]), num_cols(tdm[k]));
            }
        }
        dm.allocate_blocks();
        
        omp_for(std::size_t k, parallel::range<std::size_t>(0,dm.n_blocks()), {
            for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
                std::size_t match = half_dm[b].find_block(dm.basis().left_charge(k), dm.basis().right_charge(k));
                if (match < half_dm[b].n_blocks())
                    dm[k] += (half_dm[b][match]);
            }
        });
    }
    assert( weak_equal(dm.right_basis(), data_.right_basis()) );
    
    /// truncation
    block_matrix<Matrix, SymmGroup> U;
    block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup> S;
    truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax);
    dm = block_matrix<Matrix, SymmGroup>();
    
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, U.left_basis(), right_i, transpose(conjugate(U)), RightPaired);
    assert( mps_tensor2.reasonable() );

//    block_matrix<Matrix, SymmGroup> t = U;
//    gemm(U, S, t);
    
    block_matrix<Matrix, SymmGroup> V;
    gemm(data_, U, V);

    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, V.right_basis(), V, LeftPaired);
    assert( mps_tensor1.reasonable() );
    
    return boost::make_tuple(mps_tensor1, mps_tensor2, trunc);
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
void TwoSiteTensor<Matrix, SymmGroup>::clear()
{
    block_matrix<Matrix, SymmGroup> empty;
    swap(data(), empty);
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::swap_with(TwoSiteTensor<Matrix, SymmGroup> & b)
{
    swap(this->phys_i, b.phys_i);
    swap(this->phys_i_left, b.phys_i_left);
    swap(this->phys_i_right, b.phys_i_right);
    swap(this->left_i, b.left_i);
    swap(this->right_i, b.right_i);
    swap(this->data_, b.data_);
    swap(this->cur_storage, b.cur_storage);
    swap(this->cur_normalization, b.cur_normalization);
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> & TwoSiteTensor<Matrix, SymmGroup>::operator << (MPSTensor<Matrix, SymmGroup> const & rhs)
{
    cur_storage = TSLeftPaired;
    rhs.make_left_paired();

    // Precondition: rhs.data() and this->data() have same shape if both are left_paired
         //     assert( rhs.row_dim() == this->row_dim() &&
         // rhs.col_dim() == this->col_dim() &&
         // rhs.site_dim() == this->site_dim() );
         //     assert( rhs.data().left_basis() == this->data().left_basis() &&
         // rhs.data().right_basis() == this->data().right_basis() );
    
    left_i = rhs.row_dim();
    right_i = rhs.col_dim();
    this->data() = rhs.data();

    return *this;

}


