/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SITE_OPERATOR_H
#define SITE_OPERATOR_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"


template<class Matrix, class SymmGroup, class Dummy> class SparseOperator;

template<class Matrix, class SymmGroup>
class SiteOperator
{
    friend class SiteOperator<typename storage::constrained<Matrix>::type, SymmGroup>;
private:
    typedef typename SymmGroup::charge charge;
    typedef typename SparseOperator<Matrix, SymmGroup, void>::spin_basis_type spin_basis_type;
public:
    typedef Matrix matrix_type;
    typedef typename Matrix::size_type size_type;
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef typename maquis::traits::real_type<Matrix>::type real_type;
    typedef typename boost::ptr_vector<Matrix>::iterator block_iterator;
    typedef typename boost::ptr_vector<Matrix>::const_iterator const_block_iterator;

    SiteOperator();

    SiteOperator(Index<SymmGroup> const & rows,
                 Index<SymmGroup> const & cols);

    SiteOperator(DualIndex<SymmGroup> const & basis);

    SiteOperator(block_matrix<Matrix, SymmGroup> const&, spin_basis_type const &);

    SiteOperator& operator=(SiteOperator rhs);
    template<class OtherMatrix>
    SiteOperator& operator=(const SiteOperator<OtherMatrix, SymmGroup>& rhs);

    Index<SymmGroup> left_basis() const;
    Index<SymmGroup> right_basis() const;
    DualIndex<SymmGroup> const & basis() const;

    std::string description() const;

    Matrix &             operator[](size_type c);
    Matrix const &       operator[](size_type c) const;
    value_type &         operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c);
    value_type const &   operator()(std::pair<charge, size_type> const & r,
                                    std::pair<charge, size_type> const & c) const;
    SiteOperator &       operator+=(SiteOperator const & rhs);
    SiteOperator &       operator-=(SiteOperator const & rhs);
    SiteOperator const & operator*=(const scalar_type& v);
    SiteOperator const & operator/=(const scalar_type& v);

    size_type n_blocks() const;
    size_type find_block(charge r, charge c) const;
    bool has_block(charge r, charge c) const;
    bool has_block(std::pair<charge, size_type> const & r,
                   std::pair<charge, size_type> const & c) const;

    size_type insert_block(Matrix const &, charge, charge);
    size_type insert_block(Matrix *, charge, charge);
    void remove_block(charge r, charge c);
    void remove_block(std::size_t which);

    mutable typename parallel::scheduler_balanced_iterative::index iter_index;
    mutable typename parallel::scheduler_size_indexed::index size_index;

    real_type norm() const;
    void transpose_inplace();
    void adjoint_inplace();
    void clear();

    void match_and_add_block(Matrix const &, charge, charge);

    void reserve(charge, charge, std::size_t, std::size_t);
    void allocate_blocks();


    void resize_block(charge r, charge c,
                      size_type new_r, size_type new_c,
                      bool pretend = false);
    void resize_block(size_type pos,
                      size_type new_r, size_type new_c,
                      bool pretend = false);

    friend void swap(SiteOperator & x, SiteOperator & y)
    {
        std::swap(x.spin_, y.spin_);
        swap(x.bm_, y.bm_);
        swap(x.sparse_op, y.sparse_op);
        std::swap(x.spin_basis, y.spin_basis);
    }

    //template <class Matrix_, class SymmGroup_>
    //friend std::ostream& operator<<(std::ostream& os, SiteOperator<Matrix_, SymmGroup_> const & m);
    template <class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(symm_traits::disable_if_su2_t<SymmGroup_, std::ostream&> os,
                                    SiteOperator<Matrix_, SymmGroup_> const & m);
    template <class Matrix_, class SymmGroup_>
    friend std::ostream& operator<<(symm_traits::enable_if_su2_t<SymmGroup_, std::ostream&> os,
                                    SiteOperator<Matrix_, SymmGroup_> const & m);

    Matrix const & operator()(charge r, charge c) const
    {
        return bm_(r, c);
    }

    Matrix & operator()(charge r, charge c)
    {
        return bm_(r, c);
    }

    std::pair<const_block_iterator,const_block_iterator> blocks() const {
        return bm_.blocks();
    }

    block_matrix<Matrix, SymmGroup> & bm()
    {
        return bm_;
    }
    block_matrix<Matrix, SymmGroup> const & bm() const
    {
        return bm_;
    }

    template <class Archive>
    inline void serialize(Archive & ar, const unsigned int version);

    void update_sparse();
    SparseOperator<Matrix, SymmGroup, void> const & get_sparse() const { return sparse_op; }

    SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > & spin() { return spin_; }
    SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > const & spin() const { return spin_; }

private:
    SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > spin_;
    spin_basis_type spin_basis;
    block_matrix<Matrix, SymmGroup> bm_;
    SparseOperator<Matrix, SymmGroup, void> sparse_op;

};


#include "dmrg/block_matrix/site_operator.hpp"

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> operator*(const typename SiteOperator<Matrix,SymmGroup>::scalar_type& v,
                                          SiteOperator<Matrix, SymmGroup> bm)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> operator*(SiteOperator<Matrix, SymmGroup> bm,
                                          const typename SiteOperator<Matrix,SymmGroup>::scalar_type& v)
{
    bm *= v;
    return bm;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> operator+(SiteOperator<Matrix,SymmGroup> b1, SiteOperator<Matrix, SymmGroup> const& b2)
{
    b1 += b2;
    return b1;
}

template<class Matrix, class SymmGroup>
SiteOperator<Matrix, SymmGroup> operator-(SiteOperator<Matrix,SymmGroup> b1, SiteOperator<Matrix, SymmGroup> const& b2)
{
    b1 -= b2;
    return b1;
}

template<class Matrix, class SymmGroup>
bool shape_equal(SiteOperator<Matrix, SymmGroup> const & a, SiteOperator<Matrix, SymmGroup> const & b)
{
    return (a.basis() == b.basis() && a.spin() == b.spin());
}

template<class Matrix, class SymmGroup>
std::size_t size_of(SiteOperator<Matrix, SymmGroup> const & m)
{
    return size_of(m);
}

#include "dmrg/block_matrix/sparse_operator.h"

#endif
