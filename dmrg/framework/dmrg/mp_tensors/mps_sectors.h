/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MPS_SECTORS_H
#define MPS_SECTORS_H


#include "dmrg/block_matrix/indexing.h"
#include "dmrg/models/chem/pg_util.h"

template<class T>
T tri_min(T a, T b, T c)
{
    return std::min(std::min(a, b),
                    std::min(a, c));
}

namespace charge_detail {

    template <class SymmGroup>
    inline bool physical(typename SymmGroup::charge c) { return true; }

    template <>
    inline bool physical<SU2U1>(SU2U1::charge c) { return c[1] >= 0; }

    template <>
    inline bool physical<SU2U1PG>(SU2U1PG::charge c) { return c[1] >= 0; }
}

template <class SymmGroup>
inline std::vector<Index<SymmGroup> > allowed_sectors(std::vector<int> const& site_type,
                                                      std::vector<Index<SymmGroup> > const& phys_dims,
                                                      typename SymmGroup::charge right_end,
                                                      std::size_t Mmax)
{
    bool finitegroup = SymmGroup::finite;

    std::size_t L = site_type.size();
    
    std::vector<typename SymmGroup::charge> maximum_charges(phys_dims.size()), minimum_charges(phys_dims.size());
    for (int type=0; type<phys_dims.size(); ++type) {
        Index<SymmGroup> physc = phys_dims[type];
        physc.sort();
        maximum_charges[type] = physc.begin()->first;
        minimum_charges[type] = physc.rbegin()->first;
        if (minimum_charges[type] > maximum_charges[type]) std::swap(maximum_charges[type], minimum_charges[type]);
    }
    
    typename SymmGroup::charge maximum_total_charge=SymmGroup::IdentityCharge, minimum_total_charge=SymmGroup::IdentityCharge;
    for (int i = 0; i < L; ++i) {
        maximum_total_charge = SymmGroup::fuse(maximum_total_charge, maximum_charges[site_type[i]]);
        minimum_total_charge = SymmGroup::fuse(minimum_total_charge, minimum_charges[site_type[i]]);
    }

    Index<SymmGroup> l_triv, r_triv;
    l_triv.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
    r_triv.insert( std::make_pair(right_end, 1) );
    
    std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
    left_allowed[0] = l_triv;
    right_allowed[L] = r_triv;
    
    typename SymmGroup::charge cmaxi=maximum_total_charge, cmini=minimum_total_charge;
    for (int i = 1; i < L+1; ++i) {
        left_allowed[i] = phys_dims[site_type[i-1]] * left_allowed[i-1];
        cmaxi = SymmGroup::fuse(cmaxi, -maximum_charges[site_type[i-1]]);
        cmini = SymmGroup::fuse(cmini, -minimum_charges[site_type[i-1]]);
        
		typename Index<SymmGroup>::iterator it = left_allowed[i].begin();
        while ( it != left_allowed[i].end() )
        {
            if (!finitegroup && SymmGroup::fuse(it->first, cmaxi) < right_end)
                it = left_allowed[i].erase(it);
            else if (!finitegroup && SymmGroup::fuse(it->first, cmini) > right_end)
                it = left_allowed[i].erase(it);
            else if (!finitegroup && !charge_detail::physical<SymmGroup>(it->first))
                it = left_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
	}
    
	//--- Print left allowed terms ---//
	//std::cout << std::endl;
	//for (int i = 0; i < L+1; ++i)
	//{
	//	std::cout << "Left allowed sectors on site: " << i << std::endl;
	//	std::cout << left_allowed[i] << std::endl;
	//}
	//std::cout << std::endl;


    cmaxi=maximum_total_charge; cmini=minimum_total_charge;
    for (int i = L-1; i >= 0; --i) {
        right_allowed[i] = adjoin(phys_dims[site_type[i]]) * right_allowed[i+1];
       
        typename Index<SymmGroup>::iterator it = right_allowed[i].begin();
        while ( it != right_allowed[i].end() )
        {
            if (!finitegroup && SymmGroup::fuse(it->first, -cmaxi) > SymmGroup::IdentityCharge)
                it = right_allowed[i].erase(it);
            else if (!finitegroup && SymmGroup::fuse(it->first, -cmini) < SymmGroup::IdentityCharge)
                it = right_allowed[i].erase(it);
            else if (!finitegroup && !charge_detail::physical<SymmGroup>(it->first))
                it = right_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
        cmaxi = SymmGroup::fuse(cmaxi, -maximum_charges[site_type[i]]);
        cmini = SymmGroup::fuse(cmini, -minimum_charges[site_type[i]]);
    }
	
	//for (int i = L; i > -1; --i)
	//{
	//	std::cout << "Right allowed sectors on site: " << i << std::endl;
	//	std::cout << right_allowed[i] << std::endl;
	//}
	//std::cout << std::endl;

    for (int i = 0; i < L+1; ++i) {
        allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
        for (typename Index<SymmGroup>::iterator it = allowed[i].begin();
            it != allowed[i].end(); ++it)
            it->second = tri_min(Mmax,
                                 left_allowed[i].size_of_block(it->first),
                                 right_allowed[i].size_of_block(it->first));

		//std::cout << "Common subset of allowed sectors on site: " << i << std::endl;
		//std::cout << allowed[i] << std::endl;
    }
    
    return allowed;
}

/*
// Specialized template class for u1lpg group used for the relativistic model
template <>
inline std::vector<Index<U1LPG> > allowed_sectors<U1LPG>(std::vector<int> const& site_type,
                                                         std::vector<Index<U1LPG> > const& phys_dims,
                                                         typename U1LPG::charge right_end,
                                                         std::size_t Mmax)
{
    bool finitegroup = U1LPG::finite;
	int unbar_counter = 0;
	int bar_counter = 0;
	int delta_MK = 1;
	typename U1LPG::charge max_charge; max_charge[0] = right_end[0]/2 + 3;
	typename U1LPG::charge min_charge; min_charge[0] = right_end[0]/2 - 2;

    std::size_t L = site_type.size();
    
    std::vector<typename U1LPG::charge> maximum_charges(phys_dims.size()), minimum_charges(phys_dims.size());
    for (int type=0; type<phys_dims.size(); ++type) {
        Index<U1LPG> physc = phys_dims[type];
        physc.sort();
        maximum_charges[type] = physc.begin()->first;
        minimum_charges[type] = physc.rbegin()->first;
        if (minimum_charges[type] > maximum_charges[type]) std::swap(maximum_charges[type], minimum_charges[type]);
    }
    
    typename U1LPG::charge maximum_total_charge=U1LPG::IdentityCharge, minimum_total_charge=U1LPG::IdentityCharge;
    typename U1LPG::charge maximum_unbarred_charge=U1LPG::IdentityCharge, maximum_barred_charge=U1LPG::IdentityCharge;
    for (int i = 0; i < L; ++i) {
        maximum_total_charge = U1LPG::fuse(maximum_total_charge, maximum_charges[site_type[i]]);
        minimum_total_charge = U1LPG::fuse(minimum_total_charge, minimum_charges[site_type[i]]);
		if (site_type[i] < L/2)
			maximum_unbarred_charge = U1LPG::fuse(maximum_unbarred_charge, maximum_charges[site_type[i]]);
		else
			maximum_barred_charge = U1LPG::fuse(maximum_barred_charge, maximum_charges[site_type[i]]);
    }

    Index<U1LPG> l_triv, r_triv;
    l_triv.insert( std::make_pair(U1LPG::IdentityCharge, 1) );
    r_triv.insert( std::make_pair(right_end, 1) );
    
    std::vector<Index<U1LPG> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
    left_allowed[0] = l_triv;
    right_allowed[L] = r_triv;
    
    typename U1LPG::charge cmaxi=maximum_total_charge, cmini=minimum_total_charge;
    typename U1LPG::charge cunbar=maximum_unbarred_charge, cbar=maximum_barred_charge;
    for (int i = 1; i < L+1; ++i) {
        left_allowed[i] = phys_dims[site_type[i-1]] * left_allowed[i-1];
        cmaxi  = U1LPG::fuse(cmaxi, -maximum_charges[site_type[i-1]]);
        cmini  = U1LPG::fuse(cmini, -minimum_charges[site_type[i-1]]);

        if (site_type[i-1] < L/2) {
			++unbar_counter;
        	cunbar = U1LPG::fuse(cunbar, -maximum_charges[site_type[i-1]]);
		} else {
			++bar_counter;
        	cbar   = U1LPG::fuse(cbar, -maximum_charges[site_type[i-1]]);
		}
		if (unbar_counter > 0) {
			if (bar_counter > 0) {
				if (delta_MK > 0) {
					--delta_MK;
					max_charge[0] = max_charge[0] + 1;
				}
			}
		}

		typename Index<U1LPG>::iterator it = left_allowed[i].begin();
        while ( it != left_allowed[i].end() )
        {
            if (U1LPG::fuse(it->first, cmaxi) < max_charge)
                it = left_allowed[i].erase(it);
            else if (U1LPG::fuse(it->first, cmini) > max_charge)
                it = left_allowed[i].erase(it);
			else if ( U1LPG::fuse(it->first, cunbar) < min_charge)
				it = left_allowed[i].erase(it);
			else if ( U1LPG::fuse(it->first, cbar) < min_charge)
				it = left_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
	}
    
	//--- Print left allowed terms ---//
	std::cout << std::endl;
	for (int i = 0; i < L+1; ++i)
	{
		std::cout << "Left allowed sectors on site: " << i << std::endl;
		std::cout << left_allowed[i] << std::endl;
	}
	std::cout << std::endl;


    cmaxi=maximum_total_charge; cmini=minimum_total_charge;
	maquis::cout << maximum_total_charge << "\t" << minimum_total_charge << "\n";

    for (int i = L-1; i >= 0; --i) {
        right_allowed[i] = adjoin(phys_dims[site_type[i]]) * right_allowed[i+1];
       
        typename Index<U1LPG>::iterator it = right_allowed[i].begin();
        while ( it != right_allowed[i].end() )
        {
            if (!finitegroup && U1LPG::fuse(it->first, -cmaxi) > U1LPG::IdentityCharge)
                it = right_allowed[i].erase(it);
            else if (!finitegroup && U1LPG::fuse(it->first, -cmini) < U1LPG::IdentityCharge)
                it = right_allowed[i].erase(it);
            else if (!finitegroup && !charge_detail::physical<U1LPG>(it->first))
                it = right_allowed[i].erase(it);
            else {
                it->second = std::min(Mmax, it->second);
                ++it;
            }
        }
        cmaxi = U1LPG::fuse(cmaxi, -maximum_charges[site_type[i]]);
        cmini = U1LPG::fuse(cmini, -minimum_charges[site_type[i]]);
    }
	
	//for (int i = L; i > -1; --i)
	//{
	//	std::cout << "Right allowed sectors on site: " << i << std::endl;
	//	std::cout << right_allowed[i] << std::endl;
	//}
	//std::cout << std::endl;

    for (int i = 0; i < L+1; ++i) {
        allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
        for (typename Index<U1LPG>::iterator it = allowed[i].begin();
            it != allowed[i].end(); ++it)
            it->second = tri_min(Mmax,
                                 left_allowed[i].size_of_block(it->first),
                                 right_allowed[i].size_of_block(it->first));

		//std::cout << "Common subset of allowed sectors on site: " << i << std::endl;
		//std::cout << allowed[i] << std::endl;
    }
    return allowed;
}
*/

#endif
