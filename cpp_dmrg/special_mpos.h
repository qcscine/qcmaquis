#ifndef SPECIAL_MPOS_H
#define SPECIAL_MPOS_H

#include "mpotensor.h"

template<class Matrix>
MPOTensor<Matrix, NullGroup> identity_mpo(Index<NullGroup> phys_i)
{
    typedef Index<NullGroup>::basis_iterator bit;
    
    Index<NullGroup> aux;
    aux.insert(std::make_pair(NullGroup::Plus, 1));
    MPOTensor<Matrix, NullGroup> mpo(phys_i, aux, aux);
    mpo.multiply_by_scalar(0);
    
    for (bit up = phys_i.basis_begin(); !up.end(); ++up)
        mpo(*mpo.row_dim().basis_begin(),
            *mpo.col_dim().basis_begin(),
            *up, *up) = 1;
    
    return mpo;
}

template<class Matrix>
MPOTensor<Matrix, NullGroup> s12_sz_mpo(Index<NullGroup> phys_i)
{
    typedef Index<NullGroup>::basis_iterator bit;
    
    Index<NullGroup> aux;
    aux.insert(std::make_pair(NullGroup::Plus, 1));
    MPOTensor<Matrix, NullGroup> mpo(phys_i, aux, aux);
    mpo.multiply_by_scalar(0);
    
    int s = 1;
    for (bit up = phys_i.basis_begin(); !up.end(); ++up, s *= -1)
        mpo(*mpo.row_dim().basis_begin(),
            *mpo.col_dim().basis_begin(),
            *up, *up) = s;
    
    return mpo;
}

#endif
