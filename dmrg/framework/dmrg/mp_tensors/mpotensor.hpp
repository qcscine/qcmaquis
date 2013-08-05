/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(index_type ld,
                                        index_type rd,
                                        prempo_t const & tags,
                                        op_table_ptr tbl_)
: left_i(ld)
, right_i(rd)
, operator_table(tbl_)
, row_tags(ld, rd)
, col_tags(ld, rd)
{
    using namespace boost::tuples;
    typedef boost::tuple<index_type, index_type, tag_type, value_type> prempo_descriptor;
    typedef std::vector<prempo_descriptor> converted_prempo_t;

    if (tags.size() > 0 && operator_table.get() != NULL) {
        converted_prempo_t tmp_tags;
        
        // copy (due to const &) and convert to index_type
        for (typename prempo_t::const_iterator it = tags.begin(); it != tags.end(); ++it) {
            index_type row_i = (left_i == 1) ? 0 : index_type(get<0>(*it));
            index_type col_i = (right_i == 1) ? 0 : index_type(get<1>(*it));
            tmp_tags.push_back( prempo_descriptor(row_i, col_i, get<2>(*it), get<3>(*it)) );
        }

        std::sort(tmp_tags.begin(), tmp_tags.end(), MPOTensor_detail::row_cmp<prempo_descriptor>());

        for (typename converted_prempo_t::const_iterator it = tmp_tags.begin(); it != tmp_tags.end(); ++it)
            row_tags(get<0>(*it), get<1>(*it)) = std::make_pair(get<2>(*it), get<3>(*it));

        std::sort(tmp_tags.begin(), tmp_tags.end(), MPOTensor_detail::col_cmp<prempo_descriptor>());

        for (typename converted_prempo_t::const_iterator it = tmp_tags.begin(); it != tmp_tags.end(); ++it)
            col_tags(get<0>(*it), get<1>(*it)) = std::make_pair(get<2>(*it), get<3>(*it));
    }
    else {
        // Initialize a private operator table
        operator_table = op_table_ptr(new OPTable<Matrix, SymmGroup>());
    }
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index) const
{
    throw std::runtime_error("operator() doesn't work for MPOTensors anymore!\n");
    assert( left_index < left_i );
    assert( right_index < right_i );
    return (*operator_table)[row_tags(left_index, right_index).first];
}


template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index)
{
    throw std::runtime_error("operator() doesn't work for MPOTensors anymore!\n");
    assert( left_index < left_i );
    assert( right_index < right_i );
    typename CSRMatrix::value_type const & p = row_tags(left_index, right_index);
    return (*operator_table)[p.first];
}

template<class Matrix, class SymmGroup>
bool MPOTensor<Matrix, SymmGroup>::has(index_type left_index,
                                       index_type right_index) const
{
    assert( (row_tags.find_element(left_index, right_index) != NULL)
            == (col_tags.find_element(left_index, right_index) != NULL) );
    return row_tags.find_element(left_index, right_index) != NULL;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::set(index_type li, index_type ri, op_t const & op, value_type scale_){
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

template<class Matrix, class SymmGroup>
MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup>
MPOTensor<Matrix, SymmGroup>::at(index_type left_index, index_type right_index) const {
    typename CSRMatrix::value_type const & p = row_tags(left_index, right_index);
    return MPOTensor_detail::make_const_term_descriptor((*operator_table)[p.first], p.second);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::row_proxy MPOTensor<Matrix, SymmGroup>::row(index_type row_i) const
{  
    return row_proxy(row_tags, row_i);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::col_proxy MPOTensor<Matrix, SymmGroup>::column(index_type col_i) const
{  
    return col_proxy(col_tags, col_i);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::tag_type MPOTensor<Matrix, SymmGroup>::tag_number(index_type left_index, index_type right_index) const {
    return row_tags(left_index, right_index).first;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(const scalar_type& v)
{
    for (typename CSRMatrix::iterator1 it1 = row_tags.begin1(); it1 != row_tags.end1(); ++it1)
        for (typename CSRMatrix::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2)
            it2->second *= v; 

    for (typename CSCMatrix::iterator2 it2 = col_tags.begin2(); it2 != col_tags.end2(); ++it2)
        for (typename CSCMatrix::iterator1 it1 = it2.begin(); it1 != it2.end(); ++it1)
            it1->second *= v; 
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::divide_by_scalar(const scalar_type& v)
{
    for (typename CSRMatrix::iterator1 it1 = row_tags.begin1(); it1 != row_tags.end1(); ++it1)
        for (typename CSRMatrix::iterator2 it2 = it1.begin(); it2 != it1.end(); ++it2)
            it2->second /= v; 

    for (typename CSCMatrix::iterator2 it2 = col_tags.begin2(); it2 != col_tags.end2(); ++it2)
        for (typename CSCMatrix::iterator1 it1 = it2.begin(); it1 != it2.end(); ++it1)
            it1->second /= v; 
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::index_type MPOTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::index_type MPOTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}
