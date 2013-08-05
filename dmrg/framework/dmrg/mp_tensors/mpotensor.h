/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPOTENSOR_H
#define MPOTENSOR_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/indexing.h"
#include "utils/function_objects.h"
#include "dmrg/models/tag_table.h"

#include <iostream>
#include <set>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>

template<class Matrix, class SymmGroup>
class MPOTensor;

namespace MPOTensor_detail
{
    using namespace boost::tuples;

    template<class Matrix, class SymmGroup>
    struct pair_cmp
    {
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        bool operator()(std::pair<index_type, index_type> const & i,
                        std::pair<index_type, index_type> const & j) const
        {
            if (i.first < j.first)
                return true;
            else if (i.first > j.first)
                return false;
            else
                return i.second < j.second;
        }
    };

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

}

template <class Matrix, class SymmGroup> class column_iterator;
template <class Matrix, class SymmGroup> class compressor;
template <class Matrix, class SymmGroup> class MPOIndexer;

template<class Matrix, class SymmGroup>
class MPOTensor
{
public:
    typedef std::size_t index_type;
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;

    typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
    typedef boost::shared_ptr<OPTable<Matrix, SymmGroup> > op_table_ptr;

private:
    typedef std::pair<tag_type, value_type> internal_value_type;

public:
    typedef boost::numeric::ublas::compressed_matrix< internal_value_type,
                                                      boost::numeric::ublas::row_major
                                                      , 0, boost::numeric::ublas::unbounded_array<index_type> 
                                                    > CSRMatrix;
    typedef boost::numeric::ublas::compressed_matrix< internal_value_type,
                                                      boost::numeric::ublas::column_major
                                                      , 0, boost::numeric::ublas::unbounded_array<index_type> 
                                                    > CSCMatrix;

    typedef boost::numeric::ublas::matrix_row<const CSRMatrix> row_proxy;
    typedef boost::numeric::ublas::matrix_column<const CSCMatrix> col_proxy;

private:
    typedef std::pair<index_type, index_type> key_t;
    typedef block_matrix<Matrix, SymmGroup> value_t;
    typedef std::map<key_t, value_t, MPOTensor_detail::pair_cmp<Matrix, SymmGroup> > data_t;

    typedef std::vector<boost::tuple<std::size_t, std::size_t, tag_type, value_type> > prempo_t;
    
public:
    
    MPOTensor(index_type = 1, index_type = 1, prempo_t const & = prempo_t(), op_table_ptr = op_table_ptr());
    
    index_type row_dim() const;
    index_type col_dim() const;
    
    block_matrix<Matrix, SymmGroup> const & operator()(index_type left_index,
                                                       index_type right_index) const;
    block_matrix<Matrix, SymmGroup> & operator()(index_type left_index,
                                                 index_type right_index);

    /**** new CSR code **********/

    // tagged operator ()
    void set(index_type li, index_type ri, op_t const & op, value_type scale_ = 1.0){
        if (this->has(li, ri)) {
            assert(row_tags.find_element(li, ri)->first == col_tags.find_element(li, ri)->first);
            row_tags.find_element(li, ri)->second = scale_;
            col_tags.find_element(li, ri)->second = scale_;
            (*operator_table)[row_tags.find_element(li, ri)->first] = op;
        }
        else {
            tag_type new_tag = operator_table->register_op(op);
            row_tags(li, ri) = internal_value_type(new_tag, scale_);
            col_tags(li, ri) = internal_value_type(new_tag, scale_);
        }
    }

    // tagged operator() const
    MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup>
    at(index_type left_index, index_type right_index) const {
        typename CSRMatrix::value_type const & p = row_tags(left_index, right_index);
        return MPOTensor_detail::make_const_term_descriptor((*operator_table)[p.first], p.second);
    }

    row_proxy row(index_type row_i) const
    {
        return row_proxy(row_tags, row_i);
    }

    col_proxy column(index_type col_i) const
    {
        return col_proxy(col_tags, col_i);
    }

    tag_type tag_number(index_type left_index, index_type right_index) const {
        return row_tags(left_index, right_index).first;
    }

    bool tag_ready() const { return (row_tags.size1() > 0 && col_tags.size1() > 0 && operator_table.get() != NULL); }

    op_table_ptr get_operator_table() const { return operator_table; }

    /************************************/
    
    const data_t& data() const;

    void multiply_by_scalar(const scalar_type&);
    void divide_by_scalar(const scalar_type&);
    
    bool has(index_type left_index, index_type right_index) const;

    friend class column_iterator<Matrix, SymmGroup>;
    friend class compressor<Matrix, SymmGroup>;
    friend class MPOIndexer<Matrix, SymmGroup>;
    
private:
    data_t data_;

    CSRMatrix row_tags;
    CSCMatrix col_tags;
    op_table_ptr operator_table;
    
    index_type left_i, right_i;
};

#include "dmrg/mp_tensors/mpotensor.hpp"


#endif
