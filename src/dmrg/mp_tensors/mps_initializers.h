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
        std::size_t L = mps.length();
        
        Index<SymmGroup> l_triv, r_triv;
        l_triv.insert( std::make_pair(SymmGroup::SingletCharge, 1) );
        r_triv.insert( std::make_pair(right_end, 1) );
        
        std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
        left_allowed[0] = l_triv;
        right_allowed[L] = r_triv;
        
        for (int i = 1; i < L+1; ++i) {
            left_allowed[i] = phys * left_allowed[i-1];
            for (typename Index<SymmGroup>::iterator it = left_allowed[i].begin();
                 it != left_allowed[i].end(); ++it)
                it->second = std::min(Mmax, it->second);
        }
        for (int i = L-1; i >= 0; --i) {
            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
            for (typename Index<SymmGroup>::iterator it = right_allowed[i].begin();
                 it != right_allowed[i].end(); ++it)
                it->second = std::min(Mmax, it->second);
        }
        
        for (int i = 0; i < L+1; ++i) {
            allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
            for (typename Index<SymmGroup>::iterator it = allowed[i].begin();
                 it != allowed[i].end(); ++it)
                it->second = tri_min(Mmax,
                                     left_allowed[i].size_of_block(it->first),
                                     right_allowed[i].size_of_block(it->first));
        }
        
        for (int i = 0; i < L; ++i)
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys, allowed[i], allowed[i+1]);
        
        cout << mps.description() << endl;
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
        std::size_t L = mps.length();
        
        Index<SymmGroup> l_triv, r_triv;
        l_triv.insert( std::make_pair(SymmGroup::SingletCharge, 1) );
        r_triv.insert( std::make_pair(right_end, 1) );
        
        std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
        left_allowed[0] = l_triv;
        right_allowed[L] = r_triv;
        
        for (int i = 1; i < L+1; ++i)
            left_allowed[i] = phys * left_allowed[i-1];
        for (int i = L-1; i >= 0; --i)
            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
        
        for (int i = 0; i < L+1; ++i) {
            Index<SymmGroup> full = common_subset(left_allowed[i], right_allowed[i]);
            
            std::vector<std::size_t> sizes;
            for (typename Index<SymmGroup>::iterator it = full.begin();
                 it != full.end(); ++it) {
                sizes.push_back(std::min(left_allowed[i].size_of_block(it->first),
                                         right_allowed[i].size_of_block(it->first)));
            }
            
            std::sort(sizes.begin(), sizes.end());
            std::reverse(sizes.begin(), sizes.end());
            std::copy(sizes.begin(), sizes.end(), std::ostream_iterator<std::size_t>(cout, " ")); cout << endl;
            std::size_t threshold = sizes[std::min(sizes.size()-1, std::size_t(5))];
//            std::size_t threshold = sizes[0];
            
            for (typename Index<SymmGroup>::iterator it = full.begin();
                 it != full.end(); ++it) {
                if (std::min(left_allowed[i].size_of_block(it->first),
                             right_allowed[i].size_of_block(it->first)) >= threshold)
                    allowed[i].insert(std::make_pair(it->first, tri_min(Mmax,
                                                                        left_allowed[i].size_of_block(it->first),
                                                                        right_allowed[i].size_of_block(it->first))));
            }
        }
        
        for (int i = 0; i < L; ++i)
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys, allowed[i], allowed[i+1]);
        
        cout << mps.description() << endl;
    }
};

template<class Matrix, class SymmGroup>
struct mott_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        std::size_t L = mps.length();
        
        Index<SymmGroup> l_triv, r_triv;
        l_triv.insert( std::make_pair(SymmGroup::SingletCharge, 1) );
        r_triv.insert( std::make_pair(right_end, 1) );
        
        std::vector<Index<SymmGroup> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
        left_allowed[0] = l_triv;
        right_allowed[L] = r_triv;
        
        for (int i = 1; i < L+1; ++i) {
            left_allowed[i] = phys * left_allowed[i-1];
            for (typename Index<SymmGroup>::iterator it = left_allowed[i].begin();
                 it != left_allowed[i].end(); ++it)
                it->second = std::min(Mmax, it->second);
        }
        for (int i = L-1; i >= 0; --i) {
            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
            for (typename Index<SymmGroup>::iterator it = right_allowed[i].begin();
                 it != right_allowed[i].end(); ++it)
                it->second = std::min(Mmax, it->second);
        }
        
        for (int i = 0; i < L+1; ++i) {
            allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
            for (typename Index<SymmGroup>::iterator it = allowed[i].begin();
                 it != allowed[i].end(); ++it)
                it->second = tri_min(Mmax,
                                     left_allowed[i].size_of_block(it->first),
                                     right_allowed[i].size_of_block(it->first));
        }
        
        for (int i = 0; i < L; ++i)
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys, allowed[i], allowed[i+1], true);
        
        cout << mps.description() << endl;
    }
};

#endif
