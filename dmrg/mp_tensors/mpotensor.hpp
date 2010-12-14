#include "mp_tensors/mpotensor.h"

#include "block_matrix/reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(Index<SymmGroup> const & pd,
                                        std::size_t ld,
                                        std::size_t rd)
: data_(ld * rd, block_matrix<Matrix, SymmGroup>(pd, pd))
, left_i(ld)
, right_i(rd)
, phys_i(pd)
{ }

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_[left_index * right_i + right_index](ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type const & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index) const
{
    return data_[left_index * right_i + right_index](ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup> MPOTensor<Matrix, SymmGroup>::get_reflected() const
{
    MPOTensor<Matrix, SymmGroup> t(phys_i, right_i, left_i);
    for (std::size_t r = 0; r < left_i; ++r)
        for (std::size_t c = 0; c < right_i; ++c)
            t.data_[c*left_i+r] = data_[r*right_i+c];
    return t;
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
