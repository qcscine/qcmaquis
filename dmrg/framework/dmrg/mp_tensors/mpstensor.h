/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPSTENSOR_H
#define MPSTENSOR_H

#include <iostream>
#include <algorithm>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
//#include "solver.h"

enum boundary_flag_t {no_boundary_f,l_boundary_f,r_boundary_f};
enum MPSStorageLayout { LeftPaired, RightPaired };
// these are actually used in several places
enum Indicator { Unorm, Lnorm, Rnorm };
enum DecompMethod {QR, SVD};

static DecompMethod DefaultSolver() {return QR;} // QR or SVD

template<class Matrix, class SymmGroup>
class TwoSiteTensor;

template<class Matrix, class SymmGroup>
class MPSTensor
{
public:
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename maquis::traits::real_type<Matrix>::type real_type;
    typedef typename Matrix::value_type value_type;
    typedef double magnitude_type; // should become future (todo: Matthias, 30.04.12 / scalar-value types)
    typedef std::size_t size_type;
    using BlockMatrixType = block_matrix<Matrix, SymmGroup>;
    using BlockMatrixDiagonalType = block_matrix<typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type, SymmGroup>;

    MPSTensor(Index<SymmGroup> const & sd = Index<SymmGroup>(),
              Index<SymmGroup> const & ld = Index<SymmGroup>(),
              Index<SymmGroup> const & rd = Index<SymmGroup>(),
              bool fillrand = true,
              value_type val = 0);

    MPSTensor(Index<SymmGroup> const& sd,
              Index<SymmGroup> const& ld,
              Index<SymmGroup> const& rd,
              block_matrix<Matrix, SymmGroup> const& block,
              MPSStorageLayout layout = LeftPaired,
              Indicator = Unorm);

    Index<SymmGroup> const & site_dim() const;
    Index<SymmGroup> const & row_dim() const;
    Index<SymmGroup> const & col_dim() const;
    bool isobccompatible(Indicator) const;
    std::size_t num_elements() const;

    void replace_right_paired(block_matrix<Matrix, SymmGroup> const &, Indicator =Unorm);
    void replace_left_paired(block_matrix<Matrix, SymmGroup> const &, Indicator =Unorm);

    void add_block_to_row(block_matrix<Matrix, SymmGroup> & bm);
    void add_block_to_column(block_matrix<Matrix, SymmGroup> & bm);

    // these are not const because after a numerical test
    // they may update the status
    bool isleftnormalized(bool test = false) const;
    bool isrightnormalized(bool test = false) const;
    bool isnormalized(bool test = false) const;

    block_matrix<Matrix, SymmGroup> leftNormalizeAndReturn(DecompMethod method = DefaultSolver());
    block_matrix<Matrix, SymmGroup> rightNormalizeAndReturn(DecompMethod method = DefaultSolver());

    void leftNormalize(DecompMethod method = DefaultSolver());
    void rightNormalize(DecompMethod method = DefaultSolver());

    void shift_aux_charges(typename SymmGroup::charge);

    template <class OtherMatrix>
    void multiply_from_left(block_matrix<OtherMatrix, SymmGroup> const &);
    template <class OtherMatrix>
    void multiply_from_right(block_matrix<OtherMatrix, SymmGroup> const &);
    void multiply_by_scalar(const scalar_type&);
    void divide_by_scalar(const scalar_type&);

    scalar_type scalar_overlap(MPSTensor const &) const;
    real_type scalar_norm() const;

    // this is completely useless in C++, only exists for consistency with Python
    MPSTensor copy() const;

    BlockMatrixType & data();
    BlockMatrixType const & data() const;
    BlockMatrixType const & const_data() const;

    std::vector<block_matrix<Matrix, SymmGroup> > to_list() const;

    template<class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(std::ostream&, MPSTensor<Matrix_, SymmGroup_> const &);

    // math functions: these are not part of the Python code, but required by IETL
    MPSTensor const & operator*=(const scalar_type&);
    MPSTensor const & operator/=(const scalar_type&);

    MPSTensor const & operator+=(MPSTensor const &);
    MPSTensor const & operator-=(MPSTensor const &);

    void make_left_paired() const;
    void make_right_paired() const;

    void clear();
    void conjugate_inplace();
    void swap_with(MPSTensor & b);
    friend void swap(MPSTensor& a, MPSTensor& b){
        a.swap_with(b);
    }

    template<class Matrix_, class SymmGroup_>
    friend MPSTensor<Matrix_, SymmGroup_> join(MPSTensor<Matrix_, SymmGroup_> const &, MPSTensor<Matrix_, SymmGroup_> const &, boundary_flag_t);

    template<class Archive> void load(Archive & ar);
    template<class Archive> void save(Archive & ar) const;
    template <class Archive> void serialize(Archive & ar, const unsigned int version);

    void check_equal(MPSTensor<Matrix, SymmGroup> const &) const;
    bool reasonable() const;
    bool num_check() const; // checks for nan or inf

    Index<SymmGroup> phys_i, left_i, right_i;
private:

    // Linear algebra methods
    void performQR(BlockMatrixType& Q, BlockMatrixType& R);
    void performLQ(BlockMatrixType& L, BlockMatrixType& Q);
    void performSVD(BlockMatrixType& U, BlockMatrixType& V, BlockMatrixDiagonalType& S);

    // Class members
    mutable block_matrix<Matrix, SymmGroup> data_;
    mutable MPSStorageLayout cur_storage;
    Indicator cur_normalization;
};

// this is also required by IETL
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t,
                                       MPSTensor<Matrix, SymmGroup> m)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator*(MPSTensor<Matrix, SymmGroup> m,
                                       const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m *= t;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator/(MPSTensor<Matrix, SymmGroup> m,
                                       const typename MPSTensor<Matrix, SymmGroup>::scalar_type& t)
{
    m /= t;
    return m;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator-(MPSTensor<Matrix, SymmGroup> m,
                                       MPSTensor<Matrix, SymmGroup> const & m2)
{
    m -= m2;
    return m;
}
template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator+(MPSTensor<Matrix, SymmGroup> m,
                                       MPSTensor<Matrix, SymmGroup> const & m2)
{
    m += m2;
    return m;
}

template<class Matrix, class SymmGroup>
MPSTensor<Matrix, SymmGroup> operator-(MPSTensor<Matrix, SymmGroup> m)
{
    m *= typename MPSTensor<Matrix, SymmGroup>::scalar_type(-1.0);
    return m;
}


template<class Matrix, class SymmGroup>
std::size_t size_of(MPSTensor<Matrix, SymmGroup> const & m)
{
    return size_of(m.data());
}

#include "dmrg/mp_tensors/mpstensor.hpp"

#endif
