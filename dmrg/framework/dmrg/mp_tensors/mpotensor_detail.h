/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MPOTENSOR_DETAIL_H
#define MPOTENSOR_DETAIL_H

template<class Matrix, class SymmGroup>
class MPOTensor;

namespace MPOTensor_detail
{
    template <class Matrix, class SymmGroup>
    class term_descriptor {
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table_ptr;
    public:
        term_descriptor() {}
        term_descriptor(std::vector<std::pair<tag_type, value_type> > & term_descs, op_table_ptr op_tbl)
            : op(op_tbl->operator[](term_descs[0].first)), scale(term_descs[0].second) {}

        op_t & op;
        value_type & scale;
    };

    template <class Matrix, class SymmGroup>
    class const_term_descriptor {
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table_ptr;
    public:
        const_term_descriptor() {}
        const_term_descriptor(std::vector<std::pair<tag_type, value_type> > const & term_descs, op_table_ptr op_tbl)
            : op(op_tbl->operator[](term_descs[0].first)), scale(term_descs[0].second) {}

        op_t const & op;
        value_type const scale;
    };

    template <class ConstIterator>
    class IteratorWrapper
    {
        typedef ConstIterator internal_iterator;

    public:
        typedef IteratorWrapper<ConstIterator> self_type;
        typedef typename std::iterator_traits<internal_iterator>::value_type value_type;

        IteratorWrapper(internal_iterator i) : it_(i) { }

        void operator++() { ++it_; }
        void operator++(int) {it_++; }
        bool operator!=(self_type const & rhs) { return it_ != rhs.it_; }

        value_type index() const { return *it_; }
        value_type operator*() const {
            throw std::runtime_error("direct MPOTensor access via row iterators currently not implemented\n");
            return *it_;
        }
        
    private:
       internal_iterator it_; 
    };
    
    template <class ConstIterator>
    class row_proxy : public std::pair<ConstIterator, ConstIterator>
    {
        typedef ConstIterator internal_iterator;
        typedef std::pair<internal_iterator, internal_iterator> base;

    public:
        typedef IteratorWrapper<ConstIterator> const_iterator;
        row_proxy(internal_iterator b, internal_iterator e) : base(b, e) { } 

        const_iterator begin() const { return const_iterator(base::first); }
        const_iterator end() const { return const_iterator(base::second); }
    };

    using namespace boost::tuples;

    template<class Tuple>
    struct row_cmp
    {
        bool operator() (Tuple const & i, Tuple const & j) const
        {
            if (get<0>(i) < get<0>(j))
                return true;
            else if (get<0>(i) > get<0>(j))
                return false;
            else
                return get<1>(i) < get<1>(j);
        }
    };

    template<class Tuple>
    struct col_cmp
    {
        bool operator() (Tuple const & i, Tuple const & j) const
        {
            if (get<1>(i) < get<1>(j))
                return true;
            else if (get<1>(i) > get<1>(j))
                return false;
            else
                return get<0>(i) < get<0>(j);
        }
    };
}

#endif
