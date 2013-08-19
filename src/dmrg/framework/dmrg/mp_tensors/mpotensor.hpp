/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(std::size_t ld,
                                        std::size_t rd)
: left_i(ld)
, right_i(rd)
{ }

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::value_type & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_[std::make_pair(left_index,right_index)](ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::value_type const & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index) const
{
    return data_.find(std::make_pair(left_index,right_index))->second(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                                                                 std::size_t right_index) const
{
    assert( left_index < left_i );
    assert( right_index < right_i );
    typename data_t::const_iterator match = data_.find( std::make_pair(left_index, right_index) );
    if (match == data_.end())
        throw std::out_of_range("element not found in MPOTensor");
    return match->second;
}


template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                                                           std::size_t right_index)
{
    if (left_index > left_i)   left_i  = left_index+1;
    if (right_index > right_i) right_i = right_index+1;
    block_matrix<Matrix, SymmGroup> * ret;
    ret = &data_[std::make_pair(left_index, right_index)];
    return *ret;
}

template<class Matrix, class SymmGroup>
bool MPOTensor<Matrix, SymmGroup>::has(std::size_t left_index,
                                       std::size_t right_index) const
{
    return data_.count(std::make_pair(left_index, right_index)) > 0;
}
    
template<class Matrix, class SymmGroup>
const typename MPOTensor<Matrix, SymmGroup>::data_t& MPOTensor<Matrix, SymmGroup>::data() const {
    return this->data_;
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

template<class Matrix, class SymmGroup>
std::size_t MPOTensor<Matrix, SymmGroup>::row_dim() const
{
    return left_i;
}

template<class Matrix, class SymmGroup>
std::size_t MPOTensor<Matrix, SymmGroup>::col_dim() const
{
    return right_i;
}

