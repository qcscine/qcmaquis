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
    using std::exp; using std::pow;
    using boost::math::factorial;
    
    double w = 1.;
    for (int p=0; p<state.size(); ++p) {
        int n = boost::get<1>(state[p]);
        w *= pow(coeff[p], n) / factorial<double>(n);
    }
    return w;
}

template <class Matrix, class SymmGroup>
MPS<Matrix,SymmGroup> coherent_init(std::vector<double> const& coeff, Index<SymmGroup> const& phys,
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
        //        std::cout << "*** state:";
        //        for (int p=0; p<L; ++p)
        //            std::cout << " <" << boost::get<0>(state[p]) << " : " << boost::get<1>(state[p]) << ">";
        //        std::cout << std::endl;
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
MPS<Matrix,SymmGroup> coherent_init_dm(std::vector<double> const& coeff, Index<SymmGroup> const& phys_psi, Index<SymmGroup> const& phys_rho)
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

#endif
