/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *               2017-2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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
#include "ts_reduction.h"

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

    std::swap(data_, tmp);
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

    std::swap(data_, tmp);
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

    std::swap(data_, tmp);
    cur_storage = TSRightPaired;

    // assert( left_i == data_.left_basis() );
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> TwoSiteTensor<Matrix, SymmGroup>::make_mps() const
{
    return make_mps_(type_helper<symm_traits::HasSU2<SymmGroup>::value>());
    //return make_mps_(type_helper<false>());
}

// ALB THIS IS THE METHOD CALLED WHEN CONVERTING TWIN_MPS TO MPS TENSOR
template<class Matrix, class SymmGroup>
template<bool SU2>
MPSTensor<Matrix, SymmGroup> TwoSiteTensor<Matrix, SymmGroup>::make_mps_(type_helper<SU2>) const
{
    make_right_paired();
    return MPSTensor<Matrix, SymmGroup>(phys_i, left_i, right_i, data_, RightPaired);
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> TwoSiteTensor<Matrix, SymmGroup>::make_mps_(type_helper<true>) const
{
    make_right_paired();
    block_matrix<Matrix, SymmGroup> tmp;
    Index<SymmGroup> phys_out = ts_reduction::reduce_right(phys_i_left, phys_i_right, left_i, right_i, data_, tmp);
    return MPSTensor<Matrix, SymmGroup>(phys_out, left_i, right_i, tmp, RightPaired);
}

// ADD_NOISE_{L2R,R2L}
// ---------------------------------
// Adds noise to the DM to improve convergence
// ---------------------------------
// Note: the functions here are, unlike those in contraction::common:: specifically for TwoSiteTensors
// TODO: check whether those functions can be made to work with TwoSiteTensors too, then the two functions below may be removed

template<class Matrix, class SymmGroup>
void add_noise_l2r(block_matrix<Matrix, SymmGroup> & dm,
                MPSTensor<Matrix, SymmGroup> const & mps,
                MPOTensor<Matrix, SymmGroup> const & mpo,
                Boundary<Matrix, SymmGroup> const & boundary, // left or right depending on the direction
                double alpha)
{
    maquis::cout << "Growing, alpha = " << alpha << std::endl;

    Boundary<Matrix, SymmGroup> half_dm = contraction::Engine<Matrix, Matrix, SymmGroup>::left_boundary_tensor_mpo(mps, boundary, mpo);
    omp_for(std::size_t b, parallel::range<std::size_t>(0,half_dm.aux_dim()), {
        block_matrix<Matrix, SymmGroup> tdm;
        gemm(half_dm[b], transpose(conjugate(half_dm[b])), tdm);
        tdm *= alpha;
        std::swap(tdm, half_dm[b]);
    });

    for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
        block_matrix<Matrix, SymmGroup> const& tdm = half_dm[b];
        for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
            if (mps.data().basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
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

template<class Matrix, class SymmGroup>
void add_noise_r2l(block_matrix<Matrix, SymmGroup> & dm,
                MPSTensor<Matrix, SymmGroup> const & mps,
                MPOTensor<Matrix, SymmGroup> const & mpo,
                Boundary<Matrix, SymmGroup> const & boundary, // left or right depending on the direction
                double alpha)
{
    maquis::cout << "Growing, alpha = " << alpha << std::endl;

    Boundary<Matrix, SymmGroup> half_dm = contraction::Engine<Matrix, Matrix, SymmGroup>::right_boundary_tensor_mpo(mps, boundary, mpo);

    omp_for(std::size_t b, parallel::range<std::size_t>(0,half_dm.aux_dim()), {
        block_matrix<Matrix, SymmGroup> tdm;
        gemm(transpose(conjugate(half_dm[b])), half_dm[b], tdm);
        tdm *= alpha;
        std::swap(tdm, half_dm[b]);
    });

    for (std::size_t b = 0; b < half_dm.aux_dim(); ++b) {
        block_matrix<Matrix, SymmGroup> const& tdm = half_dm[b];
        for (std::size_t k = 0; k < tdm.n_blocks(); ++k) {
            if (mps.data().basis().has(tdm.basis().left_charge(k), tdm.basis().right_charge(k)))
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

// +-----------------+
//  SPLIT_MPS_L2R/R2L
// +-----------------+
// Methods used at the end of the optimization of a TwoSiteTensor object.
// Different for right and left sweep.
// Defines also a friend method for the case of vectors of TwoSiteTensor objects.

// -- Left sweep --

template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::split_mps_l2r(std::size_t Mmax, double cutoff, const std::vector<size_t>& keeps) const
{
    make_both_paired();

    block_matrix<Matrix, SymmGroup> u, v;
    block_diag_matrix s;

    truncation_results trunc = svd_truncate(data_, u, v, s, cutoff, Mmax, true, 0, keeps);

    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), u, LeftPaired);
    assert( mps_tensor1.reasonable() );
    gemm(s, v, u);
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, u.left_basis(), right_i, u, RightPaired);
    assert( mps_tensor2.reasonable() );
    return std::make_tuple(mps_tensor1, mps_tensor2, trunc);

}


// split_mps_l2r version which uses an externally supplied S for renormalisation
// Input: s_truncated: truncated S
//        s_full: S before truncation
//        keeps: vector of the number of eigenvalues to keep per symmetry block


template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> >
TwoSiteTensor<Matrix, SymmGroup>::split_mps_l2r(const block_diag_matrix& s_truncated,
                                                const std::vector<size_t>& keeps) const
{
    make_both_paired();

    block_matrix<Matrix, SymmGroup> u, v;

    external_truncate(data_, u, v, keeps);

    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), u, LeftPaired);
    assert( mps_tensor1.reasonable() );
    // Use s_truncated for the renormalisation
    gemm(s_truncated, v, u);
    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, u.left_basis(), right_i, u, RightPaired);
    assert( mps_tensor2.reasonable() );
    return std::make_tuple(mps_tensor1, mps_tensor2) ;
}

// -- Right sweep --

template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::split_mps_r2l(std::size_t Mmax, double cutoff, const std::vector<size_t>& keeps) const
{
    make_both_paired();

    block_matrix<Matrix, SymmGroup> u, v;
    block_diag_matrix s;

    truncation_results trunc = svd_truncate(data_, u, v, s, cutoff, Mmax, true, 0, keeps);

    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, v.left_basis(), right_i, v, RightPaired);
    assert( mps_tensor2.reasonable() );
    gemm(u, s, v);
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), v, LeftPaired);
    assert( mps_tensor1.reasonable() );
    return std::make_tuple(mps_tensor1, mps_tensor2, trunc);

}

// split_mps_r2l version which uses an externally supplied S for renormalisation
// Input: s_truncated: truncated S
//        keeps: vector of the number of eigenvalues to keep per symmetry block


template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> >
TwoSiteTensor<Matrix, SymmGroup>::split_mps_r2l(const block_diag_matrix& s_truncated, const std::vector<size_t>& keeps) const
{
    make_both_paired();

    block_matrix<Matrix, SymmGroup> u, v;

    external_truncate(data_, u, v, keeps);

    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, v.left_basis(), right_i, v, RightPaired);

    // Use s_truncated for the renormalisation
    gemm(u, s_truncated, v);
    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, u.right_basis(), v, LeftPaired);

    return std::make_tuple(mps_tensor1, mps_tensor2) ;

};

template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::predict_split_l2r(std::size_t Mmax,
                                                    double cutoff,
                                                    double alpha,
                                                    Boundary<Matrix, SymmGroup> const& left,
                                                    MPOTensor<Matrix, SymmGroup> const& mpo,
                                                    const std::vector<size_t>& keeps)
{
    make_both_paired();

    /// build reduced density matrix (with left index open)
    block_matrix<Matrix, SymmGroup> dm;
    gemm(data_, transpose(conjugate(data_)), dm, parallel::scheduler_balanced(data_));

    /// state prediction
    if (alpha != 0.) {
        Index<SymmGroup> right_phys_i = adjoin(phys_i_right) * right_i;
        MPSTensor<Matrix, SymmGroup> tmp(phys_i_left, left_i, right_phys_i, data_, LeftPaired);
        add_noise_l2r<Matrix, SymmGroup>(dm, tmp, mpo, left, alpha);
        tmp = MPSTensor<Matrix, SymmGroup>(); //why?
    }
    assert( weak_equal(dm.left_basis(), data_.left_basis()) );

    /// truncation
    block_matrix<Matrix, SymmGroup> U;
    block_diag_matrix S;
    truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax, true, keeps);
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

    return std::make_tuple(mps_tensor1, mps_tensor2, trunc);
}


template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::predict_split_r2l(std::size_t Mmax,
                                                    double cutoff,
                                                    double alpha,
                                                    Boundary<Matrix, SymmGroup> const& right,
                                                    MPOTensor<Matrix, SymmGroup> const& mpo,
                                                    const std::vector<size_t>& keeps
                                                   )
{
    make_both_paired();

    /// build reduced density matrix (with right index open)
    block_matrix<Matrix, SymmGroup> dm;
    gemm(transpose(conjugate(data_)), data_, dm, parallel::scheduler_balanced(data_));

    /// state prediction
    if (alpha != 0.) {
        Index<SymmGroup> left_phys_i = phys_i_left * left_i;
        MPSTensor<Matrix, SymmGroup> tmp(phys_i_right, left_phys_i, right_i, data_, RightPaired);
        add_noise_r2l<Matrix, SymmGroup>(dm, tmp, mpo, right, alpha);
        tmp = MPSTensor<Matrix, SymmGroup>(); //why?
    }
    assert( weak_equal(dm.right_basis(), data_.right_basis()) );

    // truncation
    block_matrix<Matrix, SymmGroup> U;
    block_diag_matrix S;
    truncation_results trunc = heev_truncate(dm, U, S, cutoff, Mmax, true, keeps);
    dm = block_matrix<Matrix, SymmGroup>();

    MPSTensor<Matrix, SymmGroup> mps_tensor2(phys_i_right, U.left_basis(), right_i, transpose(conjugate(U)), RightPaired);
    assert( mps_tensor2.reasonable() );

//    block_matrix<Matrix, SymmGroup> t = U;
//    gemm(U, S, t);

    block_matrix<Matrix, SymmGroup> V;
    gemm(data_, U, V);

    MPSTensor<Matrix, SymmGroup> mps_tensor1(phys_i_left, left_i, V.right_basis(), V, LeftPaired);
    assert( mps_tensor1.reasonable() );
    std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results> res ;
    res = std::make_tuple(mps_tensor1, mps_tensor2, trunc) ;
    return res ;
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
    std::swap(data(), empty);
}

template<class Matrix, class SymmGroup>
void TwoSiteTensor<Matrix, SymmGroup>::swap_with(TwoSiteTensor<Matrix, SymmGroup> & b)
{
    using std::swap;
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
    return operator_shift(rhs, type_helper<symm_traits::HasSU2<SymmGroup>::value>());
    //return operator_shift(rhs, type_helper<false>());
}

template<class Matrix, class SymmGroup>
template<bool SU2>
TwoSiteTensor<Matrix, SymmGroup> & TwoSiteTensor<Matrix, SymmGroup>::operator_shift(MPSTensor<Matrix, SymmGroup> const & rhs, type_helper<SU2>)
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

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> & TwoSiteTensor<Matrix, SymmGroup>::operator_shift(MPSTensor<Matrix, SymmGroup> const & rhs,
                                                                                    type_helper<true>)
{
    cur_storage = TSLeftPaired;
    rhs.make_left_paired();

    // Precondition: see above

    block_matrix<Matrix, SymmGroup> tmp;
    Index<SymmGroup> phys_out = ts_reduction::unreduce_left(phys_i_left, phys_i_right, left_i, right_i, rhs.data(), tmp);
    left_i = rhs.row_dim();
    right_i = rhs.col_dim();
    phys_i = phys_out;
    this->data() = tmp;

    return *this;
}

template<class Matrix, class SymmGroup>
std::tuple<typename TwoSiteTensor<Matrix, SymmGroup>::block_diag_matrix, truncation_results>
TwoSiteTensor<Matrix, SymmGroup>::get_S(std::size_t Mmax, double cutoff)
{
    make_both_paired();

    block_matrix<Matrix, SymmGroup> u, v;
    block_diag_matrix s;

    truncation_results trunc = svd_truncate(data_, u, v, s, cutoff, Mmax, true);

    return std::make_tuple(s, trunc);
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> const &TwoSiteTensor<Matrix, SymmGroup>::operator*=(const scalar_type& t)
{
    data() *= t;
    return *this;
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> const &TwoSiteTensor<Matrix, SymmGroup>::operator/=(const scalar_type& t)
{
    data() /= t;
    return *this;
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> const &TwoSiteTensor<Matrix, SymmGroup>::operator+=(TwoSiteTensor<Matrix, SymmGroup> const& rhs)
{
    assert( weak_equal(left_i, rhs.left_i) );
    assert( weak_equal(right_i, rhs.right_i) );
    assert( phys_i == rhs.phys_i );
    assert( phys_i_left == rhs.phys_i_left );
    assert( phys_i_right == rhs.phys_i_right );


    cur_normalization = Unorm;

    for (std::size_t i = 0; i < data().n_blocks(); ++i)
    {
        typename SymmGroup::charge lc = data().basis().left_charge(i), rc = data().basis().right_charge(i);
        if (rhs.data().has_block(lc,rc))
            data()[i] += rhs.data()(lc,rc);
    }

    return *this;
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> const &TwoSiteTensor<Matrix, SymmGroup>::operator-=(TwoSiteTensor<Matrix, SymmGroup> const& rhs)
{
    assert( weak_equal(left_i, rhs.left_i) );
    assert( weak_equal(right_i, rhs.right_i) );
    assert( phys_i == rhs.phys_i );
    assert( phys_i_left == rhs.phys_i_left );
    assert( phys_i_right == rhs.phys_i_right );


    cur_normalization = Unorm;

    for (std::size_t i = 0; i < data().n_blocks(); ++i)
    {
        typename SymmGroup::charge lc = data().basis().left_charge(i), rc = data().basis().right_charge(i);
        if (rhs.data().has_block(lc,rc))
            data()[i] -= rhs.data()(lc,rc);
    }

        return *this;
}
