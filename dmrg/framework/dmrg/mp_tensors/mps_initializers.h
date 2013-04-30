/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MPS_INITIALIZER_H
#define MPS_INITIALIZER_H

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/mp_tensors/mps_sectors.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/state_mps.h"

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
        
        maquis::cout << "Phys: " << phys << std::endl;
        maquis::cout << "Right end: " << right_end << std::endl;
        
        std::vector<Index<SymmGroup> > allowed = allowed_sectors(L, phys, right_end, Mmax);
        
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

template<class Matrix, class SymmGroup>
class coherent_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    coherent_mps_init(BaseParameters & params)
    {
        coeff = params.get<std::vector<double> >("init_coeff");
    }
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        assert(coeff.size() == mps.length());
        assert(phys.size() == 1); // only for TrivialGroup
        
        typedef typename SymmGroup::charge charge;
        
        using std::exp; using std::sqrt; using std::pow;
        using boost::math::factorial;
        
        size_t L = coeff.size();
        
        Index<SymmGroup> trivial_i;
        trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        
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
    }
private:
    std::vector<double> coeff;
};

template<class Matrix, class SymmGroup>
class coherent_dm_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    coherent_dm_mps_init(BaseParameters & params)
    {
        coeff = params.get<std::vector<double> >("init_coeff");
    }
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys_rho,
                    typename SymmGroup::charge right_end)
    {
        using std::sqrt;
        
        assert(coeff.size() == mps.length());
        assert(phys_rho.size() == 1); // only for TrivialGroup
        
        Index<SymmGroup> phys_psi;
        phys_psi.insert( std::make_pair( SymmGroup::IdentityCharge, static_cast<size_t>(sqrt(phys_rho[0].second)) ) );
        
        typedef typename SymmGroup::charge charge;
        
        using std::exp; using std::sqrt; using std::pow;
        using boost::math::factorial;
        
        size_t L = coeff.size();
        
        Index<SymmGroup> trivial_i;
        trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        
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
    }
private:
    std::vector<double> coeff;
};

template<class Matrix, class SymmGroup>
class basis_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    basis_mps_init(BaseParameters & params)
    {
        occupation = params.get<std::vector<int> >("init_basis_state");
    }
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        assert(occupation.size() == mps.length());
        assert(phys.size() == 1); // only for TrivialGroup
        typedef typename SymmGroup::charge charge;
        charge C = SymmGroup::IdentityCharge;
        
        std::vector<boost::tuple<charge, size_t> > state(mps.length());
        for (int i=0; i<mps.length(); ++i)
            state[i] = boost::make_tuple(C, occupation[i]);
        mps = state_mps<Matrix>(state, phys);
    }
    
private:
    std::vector<int> occupation;
};

template<class Matrix, class SymmGroup>
class basis_dm_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    basis_dm_mps_init(BaseParameters & params)
    {
        occupation = params.get<std::vector<int> >("init_basis_state");
    }
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys_rho,
                    typename SymmGroup::charge right_end)
    {
        assert(occupation.size() == mps.length());
        assert(phys_rho.size() == 1); // only for TrivialGroup
        typedef typename SymmGroup::charge charge;
        charge C = SymmGroup::IdentityCharge;
        
        using std::sqrt;
        size_t N = sqrt(phys_rho[0].second);
        
        std::vector<boost::tuple<charge, size_t> > state(mps.length());
        for (int i=0; i<mps.length(); ++i)
            state[i] = boost::make_tuple(C, occupation[i] + occupation[i]*N);
        mps = state_mps<Matrix>(state, phys_rho);
    }
    
private:
    std::vector<int> occupation;
};

template<class Matrix, class SymmGroup>
class basis_mps_init_generic : public mps_initializer<Matrix, SymmGroup>
{
public:
    basis_mps_init_generic(BaseParameters & params)
    {
        basis_index = params.get<std::vector<int> >("init_basis_state");
    }
    void operator()(MPS<Matrix, SymmGroup> & mps,
                    std::size_t Mmax,
                    Index<SymmGroup> const & phys,
                    typename SymmGroup::charge right_end)
    {
        assert(basis_index.size() == mps.length());
        
        std::vector<boost::tuple<typename SymmGroup::charge, size_t> > state(mps.length());
        maquis::cout << "state: ";
        for (int i=0; i<mps.length(); ++i) {
            state[i] = phys.element(basis_index[i]);
            maquis::cout << boost::get<0>(state[i]) << ":" << boost::get<1>(state[i])<< " ";
        }
        maquis::cout << "\n";
        mps = state_mps<Matrix>(state, phys);
        assert( mps[mps.length()-1].col_dim()[0].first == right_end );
    }
    
private:
    std::vector<int> basis_index;
};


#endif
