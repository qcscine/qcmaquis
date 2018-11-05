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

#ifndef TWOSITETENSOR_H
#define TWOSITETENSOR_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include <iostream>
#include <algorithm>

enum TwoSiteStorageLayout {TSRightPaired, TSLeftPaired, TSBothPaired};

// forward declaration for the declaration of friend functions to TwoSiteTensor
template<class Matrix, class SymmGroup>
class TwoSiteTensor;

template<class Matrix, class SymmGroup>
std::tuple<MPSTensor<Matrix, SymmGroup>, std::vector<MPSTensor<Matrix, SymmGroup> >, truncation_results>
    split_mps_l2r_vector(std::vector< TwoSiteTensor< Matrix, SymmGroup> > & tst_vec, std::size_t Mmax,
                        double cutoff, std::size_t Mval = 0);

template<class Matrix, class SymmGroup>
std::tuple<std::vector<MPSTensor<Matrix, SymmGroup> >, MPSTensor<Matrix, SymmGroup>, truncation_results>
    split_mps_r2l_vector(std::vector< TwoSiteTensor< Matrix, SymmGroup> > & tst_vec, std::size_t Mmax,
                        double cutoff, std::size_t Mval = 0);

// now follows the actual declaration of the TwoSiteTensor template
template<class Matrix, class SymmGroup>
class TwoSiteTensor
{
public:
    typedef std::size_t size_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename maquis::traits::real_type<Matrix>::type real_type;
    typedef typename Matrix::value_type value_type;
    typedef typename MultiIndex<SymmGroup>::index_id index_id;
    typedef typename MultiIndex<SymmGroup>::set_id set_id;
    typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type dmt;
    typedef block_matrix<dmt, SymmGroup> block_diag_matrix;

    TwoSiteTensor(MPSTensor<Matrix, SymmGroup> const & mps1,
                  MPSTensor<Matrix, SymmGroup> const & mps2);

    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    Index<SymmGroup> const & local_site_dim(short) const;

    block_matrix<Matrix, SymmGroup> & data();
    block_matrix<Matrix, SymmGroup> const & data() const;

    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, TwoSiteTensor<Matrix_, SymmGroup_> const &);

    TwoSiteTensor<Matrix, SymmGroup> & operator << ( MPSTensor<Matrix, SymmGroup> const & rhs);

    TwoSiteTensor<Matrix, SymmGroup> const & operator*=(const scalar_type&);
    TwoSiteTensor<Matrix, SymmGroup> const & operator/=(const scalar_type&);

    TwoSiteTensor<Matrix, SymmGroup> const & operator+=(TwoSiteTensor<Matrix, SymmGroup> const &);
    TwoSiteTensor<Matrix, SymmGroup> const & operator-=(TwoSiteTensor<Matrix, SymmGroup> const &);

    void make_left_paired() const;
    void make_both_paired() const;
    void make_right_paired() const;

    MPSTensor<Matrix, SymmGroup> make_mps() const;
    // get_S performs an SVD, truncates the S matrix and returns the truncated S
    std::tuple<block_diag_matrix, truncation_results> get_S(std::size_t Mmax, double cutoff);
    // -- Splitting of the TwoSiteObject with SVD --
    // split_mps_r2l & l2r versions which use an externally supplied S for renormalisation
    // Input: s_truncated: truncated S
    //        keeps: vector of the number of eigenvalues to keep per symmetry block
    std::tuple< MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results> split_mps_l2r
            (std::size_t Mmax, double cutoff, const std::vector<size_t>& keeps = std::vector<size_t>()) const;
    std::tuple< MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results> split_mps_r2l
            (std::size_t Mmax, double cutoff, const std::vector<size_t>& keeps = std::vector<size_t>()) const;
    std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> > split_mps_l2r
            (const block_diag_matrix& s_truncated, const std::vector<size_t>& keeps) const;
    std::tuple<MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup> > split_mps_r2l
            (const block_diag_matrix& s_truncated, const std::vector<size_t>& keeps) const;
    // -- Splitting of the TwoSiteObject with QR --
    std::tuple< MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results> predict_split_l2r
            (std::size_t Mmax, double cutoff, double alpha, Boundary<Matrix, SymmGroup> const& left,
             MPOTensor<Matrix, SymmGroup> const& mpo, const std::vector<size_t>& keeps = std::vector<size_t>());
    std::tuple< MPSTensor<Matrix, SymmGroup>, MPSTensor<Matrix, SymmGroup>, truncation_results> predict_split_r2l
            (std::size_t Mmax, double cutoff, double alpha, Boundary<Matrix, SymmGroup> const& right,
             MPOTensor<Matrix, SymmGroup> const& mpo, const std::vector<size_t>& keeps = std::vector<size_t>());
    // +------------------------------+
    //  Splitting of vectors of MPSTensors using SVD (sa_alg=-3)
    // +------------------------------+
    // -- Left sweep --
    friend std::tuple<MPSTensor<Matrix, SymmGroup>, std::vector<MPSTensor<Matrix, SymmGroup> >, truncation_results>
    split_mps_l2r_vector<>(std::vector<TwoSiteTensor< Matrix, SymmGroup> > & tst_vec, std::size_t Mmax,
                                double cutoff, std::size_t Mval);
    // -- Right  sweep --
    friend std::tuple< std::vector<MPSTensor<Matrix, SymmGroup> >, MPSTensor<Matrix, SymmGroup>, truncation_results>
    split_mps_r2l_vector<>(std::vector<TwoSiteTensor< Matrix, SymmGroup> > & tst_vec, std::size_t Mmax,
                         double cutoff, std::size_t Mval);

    // -- Modification methods --
    void clear();
    void swap_with(TwoSiteTensor & b);
    friend void swap(TwoSiteTensor & a, TwoSiteTensor & b) { a.swap_with(b); }
    template<class Archive> void load(Archive & ar);
    template<class Archive> void save(Archive & ar) const;
private:
    template <bool SU2> class type_helper { };
    template <bool SU2>
    MPSTensor<Matrix, SymmGroup> make_mps_(type_helper<SU2>) const;
    MPSTensor<Matrix, SymmGroup> make_mps_(type_helper<true>) const;
    template <bool SU2>
    TwoSiteTensor<Matrix, SymmGroup> & operator_shift(MPSTensor<Matrix, SymmGroup> const & rhs, type_helper<SU2>);
    TwoSiteTensor<Matrix, SymmGroup> & operator_shift(MPSTensor<Matrix, SymmGroup> const & rhs, type_helper<true>);

    MultiIndex<SymmGroup> midx;
    set_id left_paired;
    set_id right_paired;
    set_id both_paired;
    Index<SymmGroup> phys_i, phys_i_left, phys_i_right, left_i, right_i;
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable TwoSiteStorageLayout cur_storage;
    Indicator cur_normalization;
};
// this is also required by IETL
template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator*(const typename TwoSiteTensor<Matrix, SymmGroup>::scalar_type& t,
                                       TwoSiteTensor<Matrix, SymmGroup> m)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator*(TwoSiteTensor<Matrix, SymmGroup> m,
                                       const typename TwoSiteTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator/(TwoSiteTensor<Matrix, SymmGroup> m,
                                       const typename TwoSiteTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m /= t;
    return m;
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator-(TwoSiteTensor<Matrix, SymmGroup> m,
                                       TwoSiteTensor<Matrix, SymmGroup> const & m2)
{
    m -= m2;
    return m;
}
template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator+(TwoSiteTensor<Matrix, SymmGroup> m,
                                       TwoSiteTensor<Matrix, SymmGroup> const & m2)
{
    m += m2;
    return m;
}

template<class Matrix, class SymmGroup>
TwoSiteTensor<Matrix, SymmGroup> operator-(TwoSiteTensor<Matrix, SymmGroup> m)
{
    m *= typename TwoSiteTensor<Matrix, SymmGroup>::scalar_type(-1.0);
    return m;
}
#include "twositetensor.hpp"
#include "ts_ops.h"

#endif
