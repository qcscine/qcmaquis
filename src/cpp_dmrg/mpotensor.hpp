#include "mpotensor.h"

#include "reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(Index<SymmGroup> const & pd,
                                        std::size_t ld,
                                        std::size_t rd)
: data_(ld, rd, block_matrix<Matrix, SymmGroup>(pd, pd))
, phys_i(pd)
{ }

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_(left_index, right_index)(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type const & 
MPOTensor<Matrix, SymmGroup>::operator()(std::size_t left_index,
                                         std::size_t right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index) const
{
    return data_(left_index, right_index)(ket_index, bra_index);
}

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup> MPOTensor<Matrix, SymmGroup>::get_reflected() const
{
    MPOTensor<Matrix, SymmGroup> t(*this);
    t.data_ = transpose(t.data_);
    return t;
}

template<class Matrix, class SymmGroup>
void MPOTensor<Matrix, SymmGroup>::multiply_by_scalar(scalar_type v)
{
//    std::for_each(elements(data_).first, elements(data_).second, boost::lambda::_1 *= v);
    for (typename blas::general_matrix<block_matrix<Matrix, SymmGroup> >::element_iterator it = elements(data_).first;
         it != elements(data_).second; ++it)
        *it *= v;
}

template<class Matrix, class SymmGroup>
std::size_t MPOTensor<Matrix, SymmGroup>::row_dim() const
{
    return num_rows(data_);
}

template<class Matrix, class SymmGroup>
std::size_t MPOTensor<Matrix, SymmGroup>::col_dim() const
{
    return num_columns(data_);
}
