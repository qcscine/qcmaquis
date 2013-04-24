/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_COHERENT_INIT_H
#define MAQUIS_DMRG_COHERENT_INIT_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/basis_sector_iterators.h"
#include "dmrg/mp_tensors/state_mps.h"

#include <boost/tuple/tuple.hpp>

template <class SymmGroup>
double coherent_weight(std::vector<double> const& coeff, std::vector<boost::tuple<typename SymmGroup::charge, size_t> > const& state)
{
    using std::exp; using std::sqrt; using std::pow;
    using boost::math::factorial;
    
    double w = 1.;
    for (int p=0; p<state.size(); ++p) {
        int n = boost::get<1>(state[p]);
        w *= pow(coeff[p], n) * sqrt(factorial<double>(n)) / factorial<double>(n);
    }
    return w;
}

template <class Matrix, class SymmGroup>
MPS<Matrix,SymmGroup> coherent_init_join(std::vector<double> const& coeff, Index<SymmGroup> const& phys,
                                         typename SymmGroup::charge initc=SymmGroup::IdentityCharge)
{
    typedef typename SymmGroup::charge charge;
    typedef boost::tuple<charge, size_t> local_state;
    
    size_t L = coeff.size();
    
    MPS<Matrix, SymmGroup> mps;
    double prev_weight;
    bool first = true;
    basis_sector_iterator_<SymmGroup> it,end;
    for (boost::tie(it,end)=basis_sector_iterators(L, phys, initc); it!=end; ++it)
    {
        std::vector<local_state> const& state = *it;
        double weight = coherent_weight<SymmGroup>(coeff, state);
        if (mps.length() == 0) {
            mps = state_mps<Matrix>(state, phys);
        } else {
            if (first)
                mps = join(mps, state_mps<Matrix>(state, phys), prev_weight, weight);
            else
                mps = join(mps, state_mps<Matrix>(state, phys), 1., weight);
            
            first = false;
        }
        prev_weight = weight;
    }
    
    return mps;
}

template <class Matrix, class SymmGroup>
MPS<Matrix,SymmGroup> coherent_init(std::vector<double> const& coeff, Index<SymmGroup> const& phys)
{
    assert(phys.size() == 1); // only for TrivialGroup
    // TODO: require mapping phys --> dens
    
    typedef typename SymmGroup::charge charge;
    
    using std::exp; using std::sqrt; using std::pow;
    using boost::math::factorial;

    size_t L = coeff.size();
    
    Index<SymmGroup> trivial_i;
    trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));

    MPS<Matrix, SymmGroup> mps(L);
    for (int p=0; p<L; ++p) {
        int s=0;
        Matrix m(phys[s].second, 1, 0.);
        for (int ss=0; ss<phys[s].second; ++ss) {
            m(ss, 0) = pow(coeff[p], ss) * sqrt(factorial<double>(ss)) / factorial<double>(ss);
        }
        block_matrix<Matrix, SymmGroup> block;
        block.insert_block(m, SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
        
        MPSTensor<Matrix, SymmGroup> t(phys, trivial_i, trivial_i);
        t.data() = block;
        
        mps[p] = t;
    }
    return mps;
}


template <class Matrix, class SymmGroup>
MPS<Matrix,SymmGroup> coherent_init_dm_join(std::vector<double> const& coeff, Index<SymmGroup> const& phys_psi, Index<SymmGroup> const& phys_rho)
{
    typedef typename SymmGroup::charge charge;
    typedef boost::tuple<charge, size_t> local_state;
    
    size_t L = coeff.size();
    
    MPS<Matrix, SymmGroup> mps;
    double prev_weight;
    bool first = true;
    basis_sector_iterator_<SymmGroup> it1,it2,end1,end2;
    for (boost::tie(it1,end1)=basis_sector_iterators(L, phys_psi, SymmGroup::IdentityCharge); it1!=end1; ++it1)
        for (boost::tie(it2,end2)=basis_sector_iterators(L, phys_psi, SymmGroup::IdentityCharge); it2!=end2; ++it2)
    {
        std::vector<local_state> const& state1 = *it1;
        std::vector<local_state> const& state2 = *it2;
        std::vector<local_state> state_rho(L);
        
        for (int p=0; p<L; ++p) {
            boost::get<0>(state_rho[p]) = SymmGroup::IdentityCharge;
            boost::get<1>(state_rho[p]) = boost::get<1>(state1[p])*phys_psi.size_of_block(boost::get<0>(state2[p])) + boost::get<1>(state2[p]);
        }
        
        double weight = coherent_weight<SymmGroup>(coeff, state1)*coherent_weight<SymmGroup>(coeff, state2);
        if (mps.length() == 0) {
            mps = state_mps<Matrix>(state_rho, phys_rho);
        } else {
            if (first)
                mps = join(mps, state_mps<Matrix>(state_rho, phys_rho), prev_weight, weight);
            else
                mps = join(mps, state_mps<Matrix>(state_rho, phys_rho), 1., weight);
            
            first = false;
        }
        prev_weight = weight;
    }
    
    return mps;
}

template <class Matrix, class SymmGroup>
MPS<Matrix,SymmGroup> coherent_init_dm(std::vector<double> const& coeff, Index<SymmGroup> const& phys_psi, Index<SymmGroup> const& phys_rho)
{
    assert(phys_psi.size() == 1); // only for TrivialGroup
    // TODO: require mapping phys --> dens
    
    typedef typename SymmGroup::charge charge;
    
    using std::exp; using std::sqrt; using std::pow;
    using boost::math::factorial;
    
    size_t L = coeff.size();
    
    Index<SymmGroup> trivial_i;
    trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
    
    MPS<Matrix, SymmGroup> mps(L);
    for (int p=0; p<L; ++p) {
        int s=0;
        Matrix m(phys_rho[s].second, 1, 0.);
        for (int ss1=0; ss1<phys_psi[s].second; ++ss1)
            for (int ss2=0; ss2<phys_psi[s].second; ++ss2) {
            m(ss1*phys_psi[s].second+ss2, 0)  = pow(coeff[p], ss1) * sqrt(factorial<double>(ss1)) / factorial<double>(ss1);
            m(ss1*phys_psi[s].second+ss2, 0) *= pow(coeff[p], ss2) * sqrt(factorial<double>(ss2)) / factorial<double>(ss2);
        }
        block_matrix<Matrix, SymmGroup> block;
        block.insert_block(m, SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
        
        MPSTensor<Matrix, SymmGroup> t(phys_rho, trivial_i, trivial_i);
        t.data() = block;
        
        mps[p] = t;
    }
    return mps;
}


#endif
