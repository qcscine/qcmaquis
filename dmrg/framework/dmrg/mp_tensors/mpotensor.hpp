/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(index_type ld,
                                        index_type rd,
                                        prempo_t const & tags,
                                        op_table_ptr operator_table_)
: left_i(ld)
, right_i(rd)
, operator_table(operator_table_)
, row_tags(ld, rd)
, col_tags(ld, rd)
{
    using namespace boost::tuples;
    typedef std::vector<boost::tuple<index_type, index_type, tag_type, value_type> > tag_data_t;

    if (tags.size() > 0 && operator_table.get() != NULL) {
        tag_data_t tmp_tags;
        
        for (typename prempo_t::const_iterator it = tags.begin(); it != tags.end(); ++it) {

            index_type row_i = (left_i == 1) ? 0 : index_type(get<0>(*it));
            index_type col_i = (right_i == 1) ? 0 : index_type(get<1>(*it));
            boost::tuple<index_type, index_type, tag_type, value_type>
            tag = boost::make_tuple(row_i, col_i, get<2>(*it), get<3>(*it));

            tmp_tags.push_back(tag);
        }

        std::sort(tmp_tags.begin(), tmp_tags.end(), MPOTensor_detail::row_cmp<Matrix, SymmGroup>());

        for (typename tag_data_t::const_iterator it = tmp_tags.begin(); it != tmp_tags.end(); ++it)
            row_tags(get<0>(*it), get<1>(*it)) = std::make_pair(get<2>(*it), get<3>(*it));

        std::sort(tmp_tags.begin(), tmp_tags.end(), MPOTensor_detail::col_cmp<Matrix, SymmGroup>());

        for (typename tag_data_t::const_iterator it = tmp_tags.begin(); it != tmp_tags.end(); ++it)
            col_tags(get<0>(*it), get<1>(*it)) = std::make_pair(get<2>(*it), get<3>(*it));
    }
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::value_type & 
MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                         index_type right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_[std::make_pair(left_index,right_index)](ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::value_type const & 
MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                         index_type right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index) const
{
    return data_.find(std::make_pair(left_index,right_index))->second(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index) const
{
    assert( left_index < left_i );
    assert( right_index < right_i );
    return data_.at(std::make_pair(left_index, right_index));
    /*
    typename tag_data_t::const_iterator
             it = std::lower_bound(row_tags.begin(), row_tags.end(),
             boost::make_tuple(left_index, right_index, 0, 0.),
             MPOTensor_detail::row_cmp<Matrix, SymmGroup>());

    return boost::tuples::get<3>(*it) * (*operator_table)[boost::tuples::get<2>(*it)];
    */
}


template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & MPOTensor<Matrix, SymmGroup>::operator()(index_type left_index,
                                                                         index_type right_index)
{
    assert( left_index < left_i );
    assert( right_index < right_i );
    block_matrix<Matrix, SymmGroup> * ret;
    ret = &data_[std::make_pair(left_index, right_index)];
    return *ret;

    /* 
    typename tag_data_t::const_iterator
             it = std::lower_bound(row_tags.begin(), row_tags.end(),
             boost::make_tuple(left_index, right_index, 0, 0.),
             MPOTensor_detail::row_cmp<Matrix, SymmGroup>());

    return boost::tuples::get<3>(*it) * (*operator_table)[boost::tuples::get<2>(*it)];
    */
}

template<class Matrix, class SymmGroup>
bool MPOTensor<Matrix, SymmGroup>::has(index_type left_index,
                                       index_type right_index) const
{
    return data_.count(std::make_pair(left_index, right_index)) > 0;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(const scalar_type& v)
{
    for (typename std::vector<block_matrix<Matrix, SymmGroup> >::iterator it = data_.begin();
         it != data_.end(); ++it)
        *it *= v;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::divide_by_scalar(const scalar_type& v)
{
    for (typename std::vector<block_matrix<Matrix, SymmGroup> >::iterator it = data_.begin();
         it != data_.end(); ++it)
        *it /= v;
}

#ifdef AMBIENT
template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::persist() const
{
    for(typename data_t::const_iterator it = data_.begin(); it != data_.end(); ++it)
        for(index_type k = 0; k < (it->second).n_blocks(); k++)
           ambient::numeric::persist((it->second)[k]);
}
#endif

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

