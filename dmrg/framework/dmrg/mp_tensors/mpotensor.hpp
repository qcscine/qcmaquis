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
//: data_(ld * rd)
: left_i(ld)
, right_i(rd)
{ }

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_(left_index, right_index)(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type const & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         typename MPOTensor<Matrix, SymmGroup>::access_type const & bra_index) const
{
    return data_(left_index, right_index)(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> const & MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                                                                 std::size_t right_index) const
{
    assert( left_index < left_i );
    assert( right_index < right_i );
    static Timer timer("MPOTensor const lookup");
    timer.begin();
    // we need a better solution for this!
    block_matrix<Matrix, SymmGroup> * ret;
#ifdef MAQUIS_OPENMP
#pragma omp critical
#endif
    ret = &data_[std::make_pair(left_index, right_index)];
    timer.end();
    return *ret;
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> & MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                                                           std::size_t right_index)
{
    //    assert( left_index * right_i + right_index < data_.size() );
    assert( left_index < left_i );
    assert( right_index < right_i );
    static Timer timer("MPOTensor non-const lookup");
    timer.begin();
    block_matrix<Matrix, SymmGroup> * ret;
#ifdef MAQUIS_OPENMP
#pragma omp critical
#endif
//    ret = &data_[left_index][right_index];
    ret = &data_[std::make_pair(left_index, right_index)];
    timer.end();
    return *ret;
}

template<class Matrix, class SymmGroup>
bool MPOTensor<Matrix, SymmGroup>::has(std::size_t left_index,
                                       std::size_t right_index) const
{
//    static Timer timer("MPOTensor has");
//    timer.begin();
//    bool ret = data_.find(std::make_pair(left_index, right_index)) != data_.end();
//    timer.end();
//    return ret;
    return data_.count(std::make_pair(left_index, right_index)) > 0;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(scalar_type v)
{
//    std::for_each(elements(data_).first, elements(data_).second, boost::lambda::_1 *= v);
    for (typename std::vector<block_matrix<Matrix, SymmGroup> >::iterator it = data_.begin();
         it != data_.end(); ++it)
        *it *= v;
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
