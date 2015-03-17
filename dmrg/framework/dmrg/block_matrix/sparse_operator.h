/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef SPARSE_OPERATOR_H
#define SPARSE_OPERATOR_H

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"

namespace sparse_detail {

    template <class T, class SymmGroup, typename = void>
    class Entry {
    public:
        Entry();
        Entry(std::size_t r, std::size_t c, T coeff)
        : row(r), col(c), coefficient(coeff)
        {
        }

        std::size_t row, col;
        T coefficient;
    };

    template <class T, class SymmGroup>
    class Entry<T, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type> {
    public:
        typedef typename SymmGroup::subcharge subcharge;

        Entry();
        Entry(std::size_t r, std::size_t c, subcharge rspin, subcharge cspin, T coeff)
        : row(r), col(c), row_spin(rspin), col_spin(cspin), coefficient(coeff)
        {
        }

        std::size_t row, col;
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

    void update(block_matrix<Matrix, SymmGroup> const & bm, spin_basis_type const & sb)
    {
        blocks_ = std::vector<iterator>(bm.n_blocks());
        
        iterator it = data_.begin();
        for(std::size_t b = 0; b < bm.n_blocks(); ++b)
        {
            blocks_[b] = it;
            for (std::size_t ss1 = 0; ss1 < num_rows(bm[b]); ++ss1)
            for (std::size_t ss2 = 0; ss2 < num_cols(bm[b]); ++ss2)
                if (bm[b](ss1,ss2) != float_type())
                    data_.insert(it++, value_type(ss1,ss2,bm[b](ss1,ss2)));
        }
    }

private:
    std::vector<iterator> blocks_;
    std::vector<value_type> data_;
};

template<class Matrix, class SymmGroup>
class SparseOperator<Matrix, SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type> 
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
        //if(spin_basis.size() != bm.n_blocks())
        //{
        //    maquis::cout << bm << std::endl;
        //    maquis::cout << bm.n_blocks() << "  " << spin_basis.size() << std::endl;
        //    for (typename spin_basis_type::const_iterator it = spin_basis.begin(); it != spin_basis.end(); ++it)
        //    {
        //        maquis::cout << it->first.first << it->first.second << "\t" << it->second.first.size() << " " << it->second.second.size() << std::endl;
        //    }
        //    exit(1); 
        //}

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
                    //data_.push_back(value_type(ss1, ss2, left_spins[ss1], right_spins[ss2], bm[b](ss1,ss2)));
                    data_.push_back(value_type(ss1, ss2, left_spins.at(ss1), right_spins.at(ss2), bm[b](ss1,ss2)));
                    ++entry_counter;
                }
        }

        blocks_[bm.n_blocks()] = data_.size();

        //maquis::cout << "data_.size() " << data_.size() << std::endl;
        //for (int b = 0; b < blocks_.size(); ++b) {
        //    maquis::cout << "block_["<<b<<"] -> " << blocks_[b] << std::endl;
        //}
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
