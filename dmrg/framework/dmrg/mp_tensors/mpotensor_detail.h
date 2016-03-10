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

#include <boost/utility.hpp>
#include <boost/type_traits.hpp>

template<class Matrix, class SymmGroup>
class MPOTensor;

namespace MPOTensor_detail
{
    template <class T, bool C>
    struct const_type { typedef T type; };

    template <class T>
    struct const_type<T, true> { typedef const T type; };

    template <class Matrix, class SymmGroup, bool Const, typename = void>
    class term_descriptor {
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table_ptr;
    public:
        term_descriptor() {}
        term_descriptor(typename const_type<std::vector<std::pair<tag_type, value_type> >, Const>::type & term_descs,
                        op_table_ptr op_tbl)
            : op_(op_tbl->operator[](term_descs[0].first)), scale_(term_descs[0].second) {}

        typename const_type<op_t, Const>::type & op(std::size_t i=0) { return op_; }
        typename const_type<value_type, Const>::type & scale(std::size_t i=0) { return scale_; }

    private:
        typename const_type<op_t, Const>::type & op_;
        typename const_type<value_type, Const>::type & scale_;
    };

    template <class Matrix, class SymmGroup, bool Const>
    class term_descriptor<Matrix, SymmGroup, Const, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> >::type >
    {
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table_ptr;
    public:
        term_descriptor() {}
        term_descriptor(typename const_type<std::vector<std::pair<tag_type, value_type> >, Const>::type & term_descs,
                        op_table_ptr op_tbl_)
            : operator_table(op_tbl_), term_descriptors(term_descs) {}

        std::size_t size() const { return term_descriptors.size(); }
        typename const_type<op_t, Const>::type & op(std::size_t i=0) { return (*operator_table)[term_descriptors[i].first]; }
        typename const_type<value_type, Const>::type & scale(std::size_t i=0) { return term_descriptors[i].second; }
    private:
        typename const_type<std::vector<std::pair<tag_type, value_type> >, Const>::type & term_descriptors;
        op_table_ptr operator_table;
    };

    template <class ConstIterator>
    class IteratorWrapper : public std::iterator<std::forward_iterator_tag, typename std::iterator_traits<ConstIterator>::value_type>
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
