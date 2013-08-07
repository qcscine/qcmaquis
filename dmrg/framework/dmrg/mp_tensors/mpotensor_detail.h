/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename Matrix::value_type value_type;
    public:
        term_descriptor() {}
        term_descriptor(op_t & op_, value_type & s_) : op(op_), scale(s_) {}

        op_t & op;
        value_type & scale;
    };

    template <class Matrix, class SymmGroup>
    term_descriptor<Matrix, SymmGroup> make_term_descriptor(
        typename term_descriptor<Matrix, SymmGroup>::op_t & op_, typename Matrix::value_type & s_)
    {
        return term_descriptor<Matrix, SymmGroup>(op_, s_);
    }

    template <class Matrix, class SymmGroup>
    class const_term_descriptor {
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef typename Matrix::value_type value_type;
    public:
        const_term_descriptor() {}
        const_term_descriptor(op_t const & op_, value_type s_) : op(op_), scale(s_) {}

        op_t const & op;
        value_type scale;
    };

    template <class Matrix, class SymmGroup, typename Scale>
    const_term_descriptor<Matrix, SymmGroup> make_const_term_descriptor(
        block_matrix<Matrix, SymmGroup> const & op_, Scale s_)
    {
        return const_term_descriptor<Matrix, SymmGroup>(op_, s_);
    }

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
