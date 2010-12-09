#ifndef SPECIAL_MPOS_H
#define SPECIAL_MPOS_H

#include "mpotensor.h"

template<class Matrix>
MPOTensor<Matrix, NullGroup> identity_mpo(Index<NullGroup> phys_i)
{
    typedef Index<NullGroup>::basis_iterator bit;
    
    MPOTensor<Matrix, NullGroup> mpo(phys_i, 1, 1);
    mpo.multiply_by_scalar(0);
    
    for (bit up = phys_i.basis_begin(); !up.end(); ++up)
        mpo(0, 0, *up, *up) = 1;
    
    return mpo;
}

template<class Matrix>
MPOTensor<Matrix, NullGroup> s12_sz_mpo(Index<NullGroup> phys_i)
{
    typedef Index<NullGroup>::basis_iterator bit;
    
    MPOTensor<Matrix, NullGroup> mpo(phys_i, 1, 1);
    mpo.multiply_by_scalar(0);
    
    int s = 1;
    for (bit up = phys_i.basis_begin(); !up.end(); ++up, s *= -1)
        mpo(0, 0, *up, *up) = s;
    
    return mpo;
}

template<class Matrix>
std::vector<MPOTensor<Matrix, NullGroup> > s12_ising(std::size_t L, double J, double h)
{
    J *= 0.25;
    h *= 0.5;
    
    using std::make_pair;
    
    Index<NullGroup> phys; phys.insert(std::make_pair(NullGroup::Plus, 2));
    Index<NullGroup> triv; triv.insert(std::make_pair(NullGroup::Plus, 1));
    Index<NullGroup> aux; aux.insert(std::make_pair(NullGroup::Plus, 4));
    
    MPOTensor<Matrix, NullGroup> bulk(phys, 4, 4), left(phys, 1, 4), right(phys, 4, 1);
    bulk.multiply_by_scalar(0);
    left.multiply_by_scalar(0);
    right.multiply_by_scalar(0);
    
#define NGI(v) make_pair(NullGroup::Plus, v)
    
    for (int k = 0; k < 2; ++k) {
        bulk(0, 0, NGI(k), NGI(k)) = 1;
        bulk(1, 1, NGI(k), NGI(k)) = 1;
        bulk(3, 3, NGI(k), NGI(k)) = 1;
    }
        
    // Sz
    
    bulk(0, 1, NGI(0), NGI(0)) = h;
    bulk(0, 1, NGI(1), NGI(1)) = -h;
    // flip 1
    bulk(0, 2, NGI(0), NGI(1)) = J;
    bulk(0, 2, NGI(1), NGI(0)) = J;
    // flip 2
    bulk(2, 3, NGI(0), NGI(1)) = 1;
    bulk(2, 3, NGI(1), NGI(0)) = 1;
    
    left(0, 0, NGI(0), NGI(0)) = 1;
    left(0, 0, NGI(1), NGI(1)) = 1;
    left(0, 1, NGI(0), NGI(0)) = h;
    left(0, 1, NGI(1), NGI(1)) = -h;
    left(0, 2, NGI(0), NGI(1)) = J;
    left(0, 2, NGI(1), NGI(0)) = J;
    
    right(0, 0, NGI(0), NGI(0)) = h;
    right(0, 0, NGI(1), NGI(1)) = -h;
    right(2, 0, NGI(0), NGI(1)) = 1;
    right(2, 0, NGI(1), NGI(0)) = 1;
    for (int k = 0; k < 2; ++k) {
        right(1, 0, NGI(k), NGI(k)) = 1;
        right(3, 0, NGI(k), NGI(k)) = 1;
    }
    
#undef NGI
    
    std::vector<MPOTensor<Matrix, NullGroup> > ret(L, bulk);
    ret[0] = left;
    ret[L-1] = right;
    return ret;
}

#endif
