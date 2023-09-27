/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SPARSE_OPERATOR_H
#define SPARSE_OPERATOR_H

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"

namespace sparse_detail {

    template <class T, class SymmGroup, typename = void>
    class Entry {
    public:
        typedef unsigned index_type;

        Entry();
        Entry(index_type r, index_type c, T coeff)
        : row(r), col(c), coefficient(coeff)
        {
        }

        index_type row, col;
        T coefficient;
    };

    template <class T, class SymmGroup>
    class Entry<T, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>> {
    public:
        typedef typename SymmGroup::subcharge subcharge;
        typedef unsigned index_type;

        Entry();
        Entry(index_type r, index_type c, subcharge rspin, subcharge cspin, T coeff)
        : row(r), col(c), row_spin(rspin), col_spin(cspin), coefficient(coeff)
        {
        }

        index_type row, col;
        subcharge row_spin, col_spin;
        T coefficient;
    };

} // namespace sparse detail

template<class Matrix, class SymmGroup, typename = void>
class SparseOperator
{
private:
    typedef typename Matrix::value_type float_type;

public:
    typedef sparse_detail::Entry<float_type, SymmGroup> value_type;
    typedef typename std::vector<value_type>::iterator iterator;
    typedef typename std::vector<value_type>::const_iterator const_iterator;
    typedef int spin_basis_type;

    SparseOperator() {}

    SparseOperator(SiteOperator<Matrix, SymmGroup> const & bm)
    {
        update(bm);
    }

    SparseOperator(block_matrix<Matrix, SymmGroup> const & bm, spin_basis_type const & sb)
    {
        update(bm);
    }

    std::pair<const_iterator, const_iterator> block(std::size_t b) const
    {
        return std::make_pair(data_.begin() + blocks_[b], data_.begin() + blocks_[b+1]);
    }

    void update(block_matrix<Matrix, SymmGroup> const & bm, spin_basis_type const & sb)
    {
        blocks_ = std::vector<int>(bm.n_blocks() + 1);
        data_ = std::vector<value_type>();

        int entry_counter = 0;
        for(std::size_t b = 0; b < bm.n_blocks(); ++b)
        {
            blocks_[b] = entry_counter;
            for (std::size_t ss1 = 0; ss1 < num_rows(bm[b]); ++ss1)
            for (std::size_t ss2 = 0; ss2 < num_cols(bm[b]); ++ss2)
                if (bm[b](ss1,ss2) != float_type()) {
                    data_.push_back(value_type(ss1, ss2, bm[b](ss1,ss2)));
                    //data_.push_back(value_type(ss1, ss2, left_spins.at(ss1), right_spins.at(ss2), bm[b](ss1,ss2)));
                    ++entry_counter;
                }
        }

        blocks_[bm.n_blocks()] = data_.size();
    }

private:
    std::vector<int> blocks_;
    std::vector<value_type> data_;
};

template<class Matrix, class SymmGroup>
class SparseOperator<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>>
{
private:
    typedef typename Matrix::value_type float_type;
    typedef typename SymmGroup::charge charge;
    typedef std::pair<charge, charge> charge_pair;
    typedef typename SymmGroup::subcharge subcharge;

public:
    typedef sparse_detail::Entry<float_type, SymmGroup> value_type;
    typedef typename std::vector<value_type>::iterator iterator;
    typedef typename std::vector<value_type>::const_iterator const_iterator;
    typedef std::map<charge_pair, std::pair<std::vector<subcharge>, std::vector<subcharge> >, compare_pair<charge_pair> > spin_basis_type;


    SparseOperator() {}

    SparseOperator(block_matrix<Matrix, SymmGroup> const & bm, spin_basis_type const & sb)
    {
        update(bm, sb);
    }

    std::pair<const_iterator, const_iterator> block(std::size_t b) const
    {
        return std::make_pair(data_.begin() + blocks_[b], data_.begin() + blocks_[b+1]);
    }

    void update(block_matrix<Matrix, SymmGroup> const & bm, spin_basis_type const & spin_basis)
    {
        assert(spin_basis.size() >= bm.n_blocks());

        blocks_ = std::vector<int>(bm.n_blocks() + 1);
        data_ = std::vector<value_type>();

        int entry_counter = 0;
        for(std::size_t b = 0; b < bm.n_blocks(); ++b)
        {
            blocks_[b] = entry_counter;
            std::vector<subcharge> const & left_spins = spin_basis.at(std::make_pair(bm.basis().left_charge(b), bm.basis().right_charge(b))).first;
            std::vector<subcharge> const & right_spins = spin_basis.at(std::make_pair(bm.basis().left_charge(b), bm.basis().right_charge(b))).second;
            for (std::size_t ss1 = 0; ss1 < num_rows(bm[b]); ++ss1)
            for (std::size_t ss2 = 0; ss2 < num_cols(bm[b]); ++ss2)
                if (bm[b](ss1,ss2) != float_type()) {
                    data_.push_back(value_type(ss1, ss2, left_spins[ss1], right_spins[ss2], bm[b](ss1,ss2)));
                    //data_.push_back(value_type(ss1, ss2, left_spins.at(ss1), right_spins.at(ss2), bm[b](ss1,ss2)));
                    ++entry_counter;
                }
        }

        blocks_[bm.n_blocks()] = data_.size();
    }

    friend void swap(SparseOperator & x, SparseOperator & y)
    {
        swap(x.blocks_, y.blocks_);
        swap(x.data_, y.data_);
    }

private:

    std::vector<int> blocks_;
    std::vector<value_type> data_;
};

#endif
