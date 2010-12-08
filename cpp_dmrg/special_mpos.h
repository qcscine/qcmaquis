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

template<class Matrix>
std::vector<MPOTensor<Matrix, NullGroup> > s12_ising(std::size_t L, double J, double h)
{
    using std::make_pair;
    
    Index<NullGroup> phys; phys.insert(std::make_pair(NullGroup::Plus, 2));
    Index<NullGroup> triv; triv.insert(std::make_pair(NullGroup::Plus, 1));
    Index<NullGroup> aux; aux.insert(std::make_pair(NullGroup::Plus, 4));
    
    MPOTensor<Matrix, NullGroup> bulk(phys, aux, aux), left(phys, triv, aux), right(phys, aux, triv);
    bulk.multiply_by_scalar(0);
    left.multiply_by_scalar(0);
    right.multiply_by_scalar(0);
    
#define NGI(v) make_pair(NullGroup::Plus, v)
    
    for (int k = 0; k < 2; ++k) {
        bulk(NGI(0), NGI(0), NGI(k), NGI(k)) = 1;
        bulk(NGI(1), NGI(1), NGI(k), NGI(k)) = 1;
        bulk(NGI(3), NGI(3), NGI(k), NGI(k)) = 1;
    }
        
    // Sz
    bulk(NGI(0), NGI(1), NGI(0), NGI(0)) = h;
    bulk(NGI(0), NGI(1), NGI(1), NGI(1)) = -h;
    // flip 1
    bulk(NGI(0), NGI(2), NGI(0), NGI(1)) = J;
    bulk(NGI(0), NGI(2), NGI(1), NGI(0)) = J;
    // flip 2
    bulk(NGI(2), NGI(3), NGI(0), NGI(1)) = 1;
    bulk(NGI(2), NGI(3), NGI(1), NGI(0)) = 1;
    
    left(NGI(0), NGI(0), NGI(0), NGI(0)) = 1;
    left(NGI(0), NGI(0), NGI(1), NGI(1)) = 1;
    left(NGI(0), NGI(1), NGI(0), NGI(0)) = h;
    left(NGI(0), NGI(1), NGI(1), NGI(1)) = -h;
    left(NGI(0), NGI(2), NGI(0), NGI(1)) = J;
    left(NGI(0), NGI(2), NGI(1), NGI(0)) = J;
    
    right(NGI(0), NGI(0), NGI(0), NGI(0)) = h;
    right(NGI(0), NGI(0), NGI(1), NGI(1)) = -h;
    right(NGI(2), NGI(0), NGI(0), NGI(1)) = 1;
    right(NGI(2), NGI(0), NGI(1), NGI(0)) = 1;
    for (int k = 0; k < 2; ++k) {
        right(NGI(1), NGI(0), NGI(k), NGI(k)) = 1;
        right(NGI(3), NGI(0), NGI(k), NGI(k)) = 1;
    }
    
#undef NGI
    
    std::vector<MPOTensor<Matrix, NullGroup> > ret(L, bulk);
    ret[0] = left;
    ret[L-1] = right;
    return ret;
}

#endif
