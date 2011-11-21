/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SPECIAL_MPOS_H
#define SPECIAL_MPOS_H

#include "dmrg/mp_tensors/mpo.h"

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> make_mpo_identity_block(Index<SymmGroup> phys_i)
{
    block_matrix<Matrix, SymmGroup> ret;
    for (std::size_t k = 0; k < phys_i.size(); ++k)
        ret.insert_block(boost::tuples::make_tuple(Matrix(phys_i[k].second, phys_i[k].second),
                                                   phys_i[k].first,
                                                   phys_i[k].first));
    ret *= 0;
    return ret;
}

template<class Matrix, class SymmGroup>
MPOTensor<Matrix, SymmGroup> identity_mpo(Index<SymmGroup> phys_i)
{
    typedef typename Index<SymmGroup>::basis_iterator bit;
    
    MPOTensor<Matrix, SymmGroup> mpo(1, 1);
    mpo(0,0) = make_mpo_identity_block<Matrix>(phys_i);
    
    for (bit up = phys_i.basis_begin(); !up.end(); ++up)
        mpo(0, 0, *up, *up) = 1;
    
    return mpo;
}

template<class Matrix>
MPOTensor<Matrix, TrivialGroup> s12_sz_mpo(Index<TrivialGroup> phys_i)
{
    typedef Index<TrivialGroup>::basis_iterator bit;
    
    MPOTensor<Matrix, TrivialGroup> mpo(1, 1);
    mpo(0,0) = make_mpo_identity_block<Matrix>(phys_i);
    
    int s = 1;
    for (bit up = phys_i.basis_begin(); !up.end(); ++up, s *= -1)
        mpo(0, 0, *up, *up) = s;
    
    return mpo;
}

template<class Matrix>
MPO<Matrix, TrivialGroup> s12_ising(std::size_t L, double J, double h)
{
    J *= 0.25;
    h *= 0.5;
    
    using std::make_pair;
    
    Index<TrivialGroup> phys_i; phys_i.insert(std::make_pair(TrivialGroup::Plus, 2));
    Index<TrivialGroup> triv; triv.insert(std::make_pair(TrivialGroup::Plus, 1));
    Index<TrivialGroup> aux; aux.insert(std::make_pair(TrivialGroup::Plus, 4));
    
    MPOTensor<Matrix, TrivialGroup> bulk(4, 4), left(1, 4), right(4, 1);
    
    bulk(0, 0) = make_mpo_identity_block<Matrix>(phys_i);
    bulk(1, 1) = make_mpo_identity_block<Matrix>(phys_i);
    bulk(3, 3) = make_mpo_identity_block<Matrix>(phys_i);
    bulk(0, 1) = make_mpo_identity_block<Matrix>(phys_i);
    bulk(0, 2) = make_mpo_identity_block<Matrix>(phys_i);
    bulk(2, 3) = make_mpo_identity_block<Matrix>(phys_i);
    
    left(0, 0) = make_mpo_identity_block<Matrix>(phys_i);
    left(0, 1) = make_mpo_identity_block<Matrix>(phys_i);
    left(0, 2) = make_mpo_identity_block<Matrix>(phys_i);
    
    right(0, 0) = make_mpo_identity_block<Matrix>(phys_i);
    right(1, 0) = make_mpo_identity_block<Matrix>(phys_i);
    right(2, 0) = make_mpo_identity_block<Matrix>(phys_i);
    right(3, 0) = make_mpo_identity_block<Matrix>(phys_i);
    
#define NGI(v) make_pair(TrivialGroup::Plus, v)
    
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
    
    MPO<Matrix, TrivialGroup> ret(L, bulk);
    ret[0] = left;
    ret[L-1] = right;
    return ret;
}

template<class Matrix>
MPO<Matrix, U1> s12_heisenberg(std::size_t L, double Jxy, double Jz)
{
    Index<U1> phys;
    phys.insert(std::make_pair(1, 1));
    phys.insert(std::make_pair(-1, 1));
    
    block_matrix<Matrix, U1> ident, splus, sminus, sz, zero;
    
    ident.insert_block(boost::tuples::make_tuple(Matrix(1, 1, 1), -1, -1));
    ident.insert_block(boost::tuples::make_tuple(Matrix(1, 1, 1), 1, 1));
    
    splus.insert_block(boost::tuples::make_tuple(Matrix(1, 1, 1), -1, 1));
    
    sminus.insert_block(boost::tuples::make_tuple(Matrix(1, 1, 1), 1, -1));
    
    sz.insert_block(boost::tuples::make_tuple(Matrix(1, 1, 0.5), 1, 1));
    sz.insert_block(boost::tuples::make_tuple(Matrix(1, 1, -0.5), -1, -1));
    
    MPOTensor<Matrix, U1> bulk(5, 5);
    bulk(0,0) = ident;
    bulk(0,1) = Jxy/2*splus;
    bulk(1,3) = sminus;
    bulk(0,2) = Jxy/2*sminus;
    bulk(2,3) = splus;
    bulk(3,3) = ident;
    bulk(0,4) = Jz*sz;
    bulk(4,3) = sz;
    
    MPOTensor<Matrix, U1> left(1, 5);
    left(0,0) = ident;
    left(0,1) = Jxy/2*splus;
    left(0,2) = Jxy/2*sminus;
    left(0,4) = Jz*sz;
    
    MPOTensor<Matrix, U1> right(5, 1);
    right(1,0) = sminus;
    right(2,0) = splus;
    right(3,0) = ident;
    right(4,0) = sz;
    
    MPO<Matrix, U1> ret(L, bulk);
    ret[0] = left;
    ret[L-1] = right;
    
    return ret;
}

#endif
