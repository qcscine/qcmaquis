/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_INITIALIZER_H
#define MPS_INITIALIZER_H

#include "compression.h"

template<class T>
T tri_min(T a, T b, T c)
{
    return std::min(std::min(a, b),
                    std::min(a, c));
}

template<class Matrix, class SymmGroup>
struct default_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        init_sectors(mps, Mmax, phys, right_end, true);
    }
    
    void init_sectors(MPS<Matrix, SymmGroup> & mps,
                      std::size_t Mmax,
                      Index<SymmGroup> const & phys,
                      typename SymmGroup::charge right_end,
                      bool fillrand = true,
                      typename Matrix::value_type val = 0)
    {
        std::size_t L = mps.length();
        
        std::cout << "Phys: " << phys << std::endl;
        std::cout << "Right end: " << right_end << std::endl;
        
        Index<SymmGroup> physc = phys;
        physc.sort();
        typename SymmGroup::charge cmax = physc.begin()->first;
        typename SymmGroup::charge cmin = physc.rbegin()->first;
        if (cmin > cmax) std::swap(cmin, cmax);
        
        typename SymmGroup::charge cmaxL=SymmGroup::IdentityCharge, cminL=SymmGroup::IdentityCharge;
        for (int i = 1; i < L; ++i) {
            cmaxL = SymmGroup::fuse(cmaxL, cmax);
            cminL = SymmGroup::fuse(cminL, cmin);
        }
        
        Index<SymmGroup> l_triv, r_triv;
        l_triv.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
        r_triv.insert( std::make_pair(right_end, 1) );
        
        std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
        left_allowed[0] = l_triv;
        right_allowed[L] = r_triv;
        
        typename SymmGroup::charge cmaxi=cmaxL, cmini=cminL;
        for (int i = 1; i < L+1; ++i) {
            left_allowed[i] = phys * left_allowed[i-1];
            std::cout << "Bare lallowed " << i << ": " << left_allowed[i] << std::endl;
            typename Index<SymmGroup>::iterator it = left_allowed[i].begin();
            while ( it != left_allowed[i].end() )
            {
//                if (SymmGroup::fuse(it->first, cmaxi) < right_end)
//                    it = left_allowed[i].erase(it);
//                else if (SymmGroup::fuse(it->first, cmini) > right_end)
//                    it = left_allowed[i].erase(it);
//                else {
                    it->second = std::min(Mmax, it->second);
                    ++it;
//                }
            }
            cmaxi = SymmGroup::fuse(cmaxi, -cmax);
            cmini = SymmGroup::fuse(cmini, -cmin);
        }
        cmaxi=cmaxL; cmini=cminL;
        for (int i = L-1; i >= 0; --i) {
            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
            
            typename Index<SymmGroup>::iterator it = right_allowed[i].begin();
            while ( it != right_allowed[i].end() )
            {
//                if (SymmGroup::fuse(it->first, -cmaxi) > SymmGroup::IdentityCharge)
//                    it = right_allowed[i].erase(it);
//                else if (SymmGroup::fuse(it->first, -cmini) < SymmGroup::IdentityCharge)
//                    it = right_allowed[i].erase(it);
//                else {
                    it->second = std::min(Mmax, it->second);
                    ++it;
//                }
            }
            cmaxi = SymmGroup::fuse(cmaxi, -cmax);
            cmini = SymmGroup::fuse(cmini, -cmin);

        }
        
        for (int i = 0; i < L+1; ++i) {
            allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
            for (typename Index<SymmGroup>::iterator it = allowed[i].begin();
                 it != allowed[i].end(); ++it)
                it->second = tri_min(Mmax,
                                     left_allowed[i].size_of_block(it->first),
                                     right_allowed[i].size_of_block(it->first));
        }
        
        for (int i = 0; i < L; ++i) {
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys, allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        }
        
#ifndef NDEBUG
        maquis::cout << "init norm: " << norm(mps) << std::endl;
        maquis::cout << mps.description() << std::endl;
#endif
    }
};

template<class Matrix, class SymmGroup>
struct const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        default_mps_init<Matrix, SymmGroup> di;
        di.init_sectors(mps, Mmax, phys, right_end, false, 1.);
    }
};

template<class Matrix, class SymmGroup>
struct thin_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        default_mps_init<Matrix, SymmGroup> di;
        di(mps, 5, phys, right_end);
        mps = compression::l2r_compress(mps, Mmax, 1e-6); 
    }
};

template<class Matrix, class SymmGroup>
struct thin_const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        const_mps_init<Matrix, SymmGroup> di;
        di(mps, 5, phys, right_end);
        mps = compression::l2r_compress(mps, Mmax, 1e-6); 
    }
};


template<class Matrix, class SymmGroup>
struct empty_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        std::size_t L = mps.length();
        
        Index<SymmGroup> triv;
        triv.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
        
        for (int i = 0; i < L; ++i)
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys, triv, triv);        
    }
};

template<class Matrix>
struct linear_mps_init : public mps_initializer<Matrix, U1>
{
    
    void operator()(MPS<Matrix, U1> & mps,
                    std::size_t Mmax,
                    Index<U1> const & phys,
                    U1::charge right_end)
    {
        init_sectors(mps, Mmax, phys, right_end, true);
    }
    
    void init_sectors(MPS<Matrix, U1> & mps,
                      std::size_t Mmax,
                      Index<U1> const & phys,
                      U1::charge right_end,
                      bool fillrand = true,
                      typename Matrix::value_type val = 1.)
    {
        std::size_t L = mps.length();
        int delta = 1;
        
        Index<U1> physc = phys;
        physc.sort();
        U1::charge cmax = physc.begin()->first;
        U1::charge cmin = physc.rbegin()->first;
        if (cmin > cmax) std::swap(cmin, cmax);
        
        U1::charge cmaxL=U1::IdentityCharge, cminL=U1::IdentityCharge;
        for (int i = 1; i < L; ++i) {
            cmaxL = U1::fuse(cmaxL, cmax);
            cminL = U1::fuse(cminL, cmin);
        }
        
        Index<U1> l_triv, r_triv;
        l_triv.insert( std::make_pair(U1::IdentityCharge, 1) );
        r_triv.insert( std::make_pair(right_end, 1) );
        
        std::vector<Index<U1> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
        left_allowed[0] = l_triv;
        right_allowed[L] = r_triv;
        
        for (int i = 1; i < L+1; ++i) {
            U1::charge cell_end = ((i-1)*right_end) / L + 1;
            int num_site_right = L/right_end - (i-1)%(L/right_end) - 1;
            U1::charge cmaxi = num_site_right*cmax;
            U1::charge cmini = num_site_right*cmin;
            
            left_allowed[i] = phys * left_allowed[i-1];
            Index<U1>::iterator it = left_allowed[i].begin();
            while ( it != left_allowed[i].end() )
            {
                if (U1::fuse(it->first, cmaxi) +delta < cell_end)
                    it = left_allowed[i].erase(it);
                else if (U1::fuse(it->first, cmini) -delta > cell_end)
                    it = left_allowed[i].erase(it);
                else {
                    it->second = std::min(Mmax, it->second);
                    ++it;
                }
            }
        }
        
        for (int i = L-1; i >= 0; --i) {
            U1::charge cell_end = (i*right_end) / L;
            int num_site_left = i%(L/right_end);
            U1::charge cmaxi = num_site_left*(cmax);
            U1::charge cmini = num_site_left*(cmin);
            
            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
            Index<U1>::iterator it = right_allowed[i].begin();
            while ( it != right_allowed[i].end() )
            {
                if (U1::fuse(it->first, -cmaxi) -delta > cell_end)
                    it = right_allowed[i].erase(it);
                else if (U1::fuse(it->first, -cmini) +delta < cell_end)
                    it = right_allowed[i].erase(it);
                else {
                    it->second = std::min(Mmax, it->second);
                    ++it;
                }
            }
        }
        
        for (int i = 0; i < L+1; ++i) {
            allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
            for (Index<U1>::iterator it = allowed[i].begin();
                 it != allowed[i].end(); ++it)
                it->second = tri_min(Mmax,
                                     left_allowed[i].size_of_block(it->first),
                                     right_allowed[i].size_of_block(it->first));
        }
        
        for (int i = 0; i < L; ++i) {
            mps[i] = MPSTensor<Matrix, U1>(phys, allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        }
        
#ifndef NDEBUG
        maquis::cout << "init norm: " << norm(mps) << std::endl;
        maquis::cout << mps.description() << std::endl;
#endif
    }
};



#endif
