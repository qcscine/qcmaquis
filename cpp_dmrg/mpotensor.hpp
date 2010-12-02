#include "mpotensor.h"

#include "reshapes.h"

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup>::MPOTensor(Index<SymmGroup> const & sd,
                                        Index<SymmGroup> const & ld,
                                        Index<SymmGroup> const & rd)
: phys_i(sd)
, left_i(ld)
, right_i(rd)
, data_(sd*ld, sd*rd)
, cur_storage(LeftUp)
, cur_normalization(U)
{
    data_.fill_with_random(drand48);
}

template<class Matrix, class SymmGroup>
typename MPOTensor<Matrix, SymmGroup>::scalar_type & 
MPOTensor<Matrix, SymmGroup>::operator()(MPOTensor<Matrix, SymmGroup>::access_type const & left_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & right_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & ket_index,
                                         MPOTensor<Matrix, SymmGroup>::access_type const & bra_index)
{
    return data_(calculate_index<SymmGroup, 2>((phys_i, left_i),
                                               (ket_index, left_index)),
                 calculate_index<SymmGroup, 2>((phys_i, right_i),
                                               (bra_index, right_index)));
}
