/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MPS_INITIALIZER_H
#define MPS_INITIALIZER_H

#include <fstream>
#include <sstream>
#include <algorithm>

#include <boost/tokenizer.hpp>

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/mp_tensors/mps_sectors.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/state_mps.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers_helper.h"

// ========================================
//  IMPLEMENTATION OF THE MPS INITIALIZERS
// ========================================

// == DEFAULT_MPS_INIT ==
// Note that this method supports both random and constant initialization of the MPS 

template<class Matrix, class SymmGroup>
struct default_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    /** 
     * @brief Class constructor 
     */
    default_mps_init(BaseParameters & parms, std::vector<Index<SymmGroup> > const& phys_dims_,
                     typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
    : init_bond_dimension(parms["init_bond_dimension"])
    , phys_dims(phys_dims_)
    , right_end(right_end_)
    , site_type(site_type_)
    { }
    
    /**
     * @brief Functor operator called to generate the MPS
     * 
     * Note that fillrand variable is set here to true - so random initialization is done
     * by default.
     * 
     * @param mps output MPS
     */
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        init_sectors(mps, this->init_bond_dimension, true);
    }
    
    // Main routine
    void init_sectors(MPS<Matrix, SymmGroup> & mps, size_t Mmax, bool fillrand=true, typename Matrix::value_type val=0)
    {
        parallel::scheduler_balanced scheduler(mps.length());
        std::size_t L = mps.length();
        maquis::cout << "Right end: " << right_end << std::endl;
        // Compute the indexes which are allowed by symmetry
        std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, Mmax);
        // Populates the MPS tensor
        omp_for(size_t i, parallel::range<size_t>(0,L), {
            parallel::guard proc(scheduler(i));
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], fillrand, val);
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        });
        
#ifndef NDEBUG
        maquis::cout << "init norm: " << norm(mps) << std::endl;
        maquis::cout << mps.description() << std::endl;
#endif
    }

    // Class attributes
    int init_bond_dimension;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
};

/**
 * @brief Const MPS initializer
 *
 * This initializer populates the MPS with all 1s
 */
template<class Matrix, class SymmGroup>
struct const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    /** @brief Class costructor */
    const_mps_init(BaseParameters & parms,
                   std::vector<Index<SymmGroup> > const& phys_dims,
                   typename SymmGroup::charge right_end,
                   std::vector<int> const& site_type)
    : di(parms, phys_dims, right_end, site_type)
    { }

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        di.init_sectors(mps, di.init_bond_dimension, false, 1.);
    }

    default_mps_init<Matrix, SymmGroup> di;
};

template<class Matrix, class SymmGroup>
struct thin_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    thin_mps_init(BaseParameters & parms,
                  std::vector<Index<SymmGroup> > const& phys_dims,
                  typename SymmGroup::charge right_end,
                  std::vector<int> const& site_type)
    : di(parms, phys_dims, right_end, site_type)
    { }

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        di.init_sectors(mps, 5, true);
        mps = compression::l2r_compress(mps, di.init_bond_dimension, 1e-6);
    }

    default_mps_init<Matrix, SymmGroup> di;
};

template<class Matrix, class SymmGroup>
struct thin_const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    thin_const_mps_init(BaseParameters & parms,
                        std::vector<Index<SymmGroup> > const& phys_dims,
                        typename SymmGroup::charge right_end,
                        std::vector<int> const& site_type)
    : di(parms, phys_dims, right_end, site_type)
    { }

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        di.init_sectors(mps, 5, false, 1.);
        mps = compression::l2r_compress(mps, di.init_bond_dimension, 1e-6);
    }

    default_mps_init<Matrix, SymmGroup> di;
};


template<class Matrix, class SymmGroup>
struct empty_mps_init : public mps_initializer<Matrix, SymmGroup>
{
    empty_mps_init(std::vector<Index<SymmGroup> > const& phys_dims_,
                   std::vector<int> const& site_type_)
    : phys_dims(phys_dims_)
    , site_type(site_type_)
    { }
    
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        std::size_t L = mps.length();
        
        Index<SymmGroup> triv;
        triv.insert( std::make_pair(SymmGroup::IdentityCharge, 1) );
        
        for (int i = 0; i < L; ++i)
            mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], triv, triv);
    }
    
private:
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_type;
};

//template<class Matrix>
//struct linear_mps_init : public mps_initializer<Matrix, U1>
//{
//    
//    void operator()(MPS<Matrix, U1> & mps,
//                    std::size_t Mmax,
//                    Index<U1> const & phys,
//                    U1::charge right_end)
//    {
//        init_sectors(mps, Mmax, phys, right_end, true);
//    }
//    
//    void init_sectors(MPS<Matrix, U1> & mps,
//                      std::size_t Mmax,
//                      Index<U1> const & phys,
//                      U1::charge right_end,
//                      bool fillrand = true,
//                      typename Matrix::value_type val = 1.)
//    {
//        std::size_t L = mps.length();
//        int delta = 1;
//        
//        Index<U1> physc = phys;
//        physc.sort();
//        U1::charge cmax = physc.begin()->first;
//        U1::charge cmin = physc.rbegin()->first;
//        if (cmin > cmax) std::swap(cmin, cmax);
//        
//        U1::charge cmaxL=U1::IdentityCharge, cminL=U1::IdentityCharge;
//        for (int i = 1; i < L; ++i) {
//            cmaxL = U1::fuse(cmaxL, cmax);
//            cminL = U1::fuse(cminL, cmin);
//        }
//        
//        Index<U1> l_triv, r_triv;
//        l_triv.insert( std::make_pair(U1::IdentityCharge, 1) );
//        r_triv.insert( std::make_pair(right_end, 1) );
//        
//        std::vector<Index<U1> > left_allowed(L+1), right_allowed(L+1), allowed(L+1);
//        left_allowed[0] = l_triv;
//        right_allowed[L] = r_triv;
//        
//        for (int i = 1; i < L+1; ++i) {
//            U1::charge cell_end = ((i-1)*right_end) / L + 1;
//            int num_site_right = L/right_end - (i-1)%(L/right_end) - 1;
//            U1::charge cmaxi = num_site_right*cmax;
//            U1::charge cmini = num_site_right*cmin;
//            
//            left_allowed[i] = phys * left_allowed[i-1];
//            Index<U1>::iterator it = left_allowed[i].begin();
//            while ( it != left_allowed[i].end() )
//            {
//                if (U1::fuse(it->first, cmaxi) +delta < cell_end)
//                    it = left_allowed[i].erase(it);
//                else if (U1::fuse(it->first, cmini) -delta > cell_end)
//                    it = left_allowed[i].erase(it);
//                else {
//                    it->second = std::min(Mmax, it->second);
//                    ++it;
//                }
//            }
//        }
//        
//        for (int i = L-1; i >= 0; --i) {
//            U1::charge cell_end = (i*right_end) / L;
//            int num_site_left = i%(L/right_end);
//            U1::charge cmaxi = num_site_left*(cmax);
//            U1::charge cmini = num_site_left*(cmin);
//            
//            right_allowed[i] = adjoin(phys) * right_allowed[i+1];
//            Index<U1>::iterator it = right_allowed[i].begin();
//            while ( it != right_allowed[i].end() )
//            {
//                if (U1::fuse(it->first, -cmaxi) -delta > cell_end)
//                    it = right_allowed[i].erase(it);
//                else if (U1::fuse(it->first, -cmini) +delta < cell_end)
//                    it = right_allowed[i].erase(it);
//                else {
//                    it->second = std::min(Mmax, it->second);
//                    ++it;
//                }
//            }
//        }
//        
//        for (int i = 0; i < L+1; ++i) {
//            allowed[i] = common_subset(left_allowed[i], right_allowed[i]);
//            for (Index<U1>::iterator it = allowed[i].begin();
//                 it != allowed[i].end(); ++it)
//                it->second = tri_min(Mmax,
//                                     left_allowed[i].size_of_block(it->first),
//                                     right_allowed[i].size_of_block(it->first));
//        }
//        
//        for (int i = 0; i < L; ++i) {
//            mps[i] = MPSTensor<Matrix, U1>(phys, allowed[i], allowed[i+1], fillrand, val);
//            mps[i].divide_by_scalar(mps[i].scalar_norm());
//        }
//        
//#ifndef NDEBUG
//        maquis::cout << "init norm: " << norm(mps) << std::endl;
//        maquis::cout << mps.description() << std::endl;
//#endif
//    }
//};

template<class Matrix, class SymmGroup>
class coherent_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    coherent_mps_init(BaseParameters & params,
                      std::vector<Index<SymmGroup> > const& phys_dims_,
                      std::vector<int> const& site_type_)
    : coeff(params["init_coeff"].as<std::vector<double> >())
    , phys_dims(phys_dims_)
    , site_type(site_type_)
    { }
    
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        assert(coeff.size() == mps.length());
        if (phys_dims[0].size() != 1) throw std::runtime_error("coherent_mps_init only for TrivialGroup.");
        
        typedef typename SymmGroup::charge charge;
        
        using std::exp; using std::sqrt; using std::pow;
        using boost::math::factorial;
        
        size_t L = coeff.size();
        
        Index<SymmGroup> trivial_i;
        trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        
        for (int p=0; p<L; ++p) {
            int s=0;
            Matrix m(phys_dims[site_type[p]][s].second, 1, 0.);
            for (int ss=0; ss<phys_dims[site_type[p]][s].second; ++ss) {
                m(ss, 0) = pow(coeff[p], ss) * sqrt(factorial<double>(ss)) / factorial<double>(ss);
            }
            block_matrix<Matrix, SymmGroup> block;
            block.insert_block(m, SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
            
            MPSTensor<Matrix, SymmGroup> t(phys_dims[site_type[p]], trivial_i, trivial_i);
            t.data() = block;
            
            swap(mps[p], t);
        }
    }
private:
    std::vector<double> coeff;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_type;
};

template<class Matrix, class SymmGroup>
class coherent_dm_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    coherent_dm_mps_init(BaseParameters & params,
                         std::vector<Index<SymmGroup> > const& phys_dims_,
                         std::vector<int> const& site_type_)
    : coeff(params["init_coeff"].as<std::vector<double> >())
    , phys_rho_dims(phys_dims_)
    , site_type(site_type_)
    { }

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        using std::sqrt;
        
        assert(coeff.size() == mps.length());
        if (phys_rho_dims[0].size() != 1) throw std::runtime_error("coherent_dm_mps_init only for TrivialGroup.");
        
        std::vector<Index<SymmGroup> > phys_psi_dims(phys_rho_dims.size());
        for (int type=0; type<phys_rho_dims.size(); ++type)
            phys_psi_dims[type].insert( std::make_pair( SymmGroup::IdentityCharge, static_cast<size_t>(sqrt(double(phys_rho_dims[type][0].second))) ) );
        
        typedef typename SymmGroup::charge charge;
        
        using std::exp; using std::sqrt; using std::pow;
        using boost::math::factorial;
        
        size_t L = coeff.size();
        
        Index<SymmGroup> trivial_i;
        trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
        
        for (int p=0; p<L; ++p) {
            int s=0;
            Matrix m(phys_rho_dims[site_type[p]][s].second, 1, 0.);
            for (int ss1=0; ss1<phys_psi_dims[site_type[p]][s].second; ++ss1)
                for (int ss2=0; ss2<phys_psi_dims[site_type[p]][s].second; ++ss2) {
                    m(ss1*phys_psi_dims[site_type[p]][s].second+ss2, 0)  = pow(coeff[p], ss1) * sqrt(factorial<double>(ss1)) / factorial<double>(ss1);
                    m(ss1*phys_psi_dims[site_type[p]][s].second+ss2, 0) *= pow(coeff[p], ss2) * sqrt(factorial<double>(ss2)) / factorial<double>(ss2);
                }
            block_matrix<Matrix, SymmGroup> block;
            block.insert_block(m, SymmGroup::IdentityCharge, SymmGroup::IdentityCharge);
            
            MPSTensor<Matrix, SymmGroup> t(phys_rho_dims[site_type[p]], trivial_i, trivial_i);
            t.data() = block;
            
            swap(mps[p], t);
        }
    }
private:
    std::vector<double> coeff;
    std::vector<Index<SymmGroup> > phys_rho_dims;
    std::vector<int> site_type;
};

template<class Matrix, class SymmGroup>
class basis_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    basis_mps_init(BaseParameters & params,
                   std::vector<Index<SymmGroup> > const& phys_dims_,
                   std::vector<int> const& site_type_)

    : occupation(params["init_basis_state"].as<std::vector<int> >())
    , phys_dims(phys_dims_)
    , site_type(site_type_)
    { }

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        assert(occupation.size() == mps.length());
        if (phys_dims[0].size() != 1) throw std::runtime_error("basis_mps_init only for TrivialGroup.");
        typedef typename SymmGroup::charge charge;
        charge C = SymmGroup::IdentityCharge;
        
        std::vector<boost::tuple<charge, size_t> > state(mps.length());
        for (int i=0; i<mps.length(); ++i)
            state[i] = boost::make_tuple(C, occupation[i]);
        mps = state_mps<Matrix>(state, phys_dims, site_type);
    }
    
private:
    std::vector<int> occupation;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_type;
};

template<class Matrix, class SymmGroup>
class basis_dm_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    basis_dm_mps_init(BaseParameters & params,
                      std::vector<Index<SymmGroup> > const& phys_dims_,
                      std::vector<int> const& site_type_)
    : occupation(params["init_basis_state"].as<std::vector<int> >())
    , phys_rho_dims(phys_dims_)
    , site_type(site_type_)
    { }
    
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        assert(occupation.size() == mps.length());
        if (phys_rho_dims.size() != 1) throw std::runtime_error("basis_dm_mps_init only for unique site basis.");
        if (phys_rho_dims[0].size() != 1) throw std::runtime_error("basis_dm_mps_init only for TrivialGroup.");
        typedef typename SymmGroup::charge charge;
        charge C = SymmGroup::IdentityCharge;
        
        using std::sqrt;
        size_t N = sqrt(double(phys_rho_dims[0][0].second));
        
        std::vector<boost::tuple<charge, std::size_t> > state(mps.length());
        for (int i=0; i<mps.length(); ++i)
            state[i] = boost::make_tuple(C, occupation[i] + occupation[i]*N);
        mps = state_mps<Matrix>(state, phys_rho_dims, site_type);
    }
    
private:
    std::vector<int> occupation;
    std::vector<Index<SymmGroup> > phys_rho_dims;
    std::vector<int> site_type;
};

/**
 * @brief ONV MPS initializer
 *
 * The MPS is initialized from a unique ONV.
 * The initial MPS has, therefore, m=1
 */
template<class Matrix, class SymmGroup>
class basis_mps_init_generic : public mps_initializer<Matrix, SymmGroup>
{
public:
    // Types definition
    typedef std::vector<boost::tuple<typename SymmGroup::charge, size_t> > state_type;

    /** 
     * @brief Class constructor from a parameter object
     * @param params Parameter container.
     * @param phys_dims_ Vector with the physical index per site type.
     * @param right_end_ Overall symmetry sector to which the MPS belongs.
     * @param site_type_ Vector with size == the lattice size, with the type of each site.
     */
    basis_mps_init_generic(BaseParameters & params, const std::vector<Index<SymmGroup> >& phys_dims_,
                           typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
        : basis_index(params["init_basis_state"].as<std::vector<int> >()), phys_dims(phys_dims_),
          right_end(right_end_), site_type(site_type_)
    { }
    
    /** @brief Operator (), called when the MPS is constructed */
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        // assert(basis_index.size() == mps.length());
        auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(basis_index, phys_dims, site_type, mps.length());
        mps = state_mps<Matrix>(state, phys_dims, site_type);
#ifndef NDEBUG
        for (int i = 0 ; i < basis_index.size() ; i++ ) {
          maquis::cout << "state: ";
          maquis::cout << boost::get<0>(state[i]) << ":" << boost::get<1>(state[i])<< " ";
          maquis::cout << "\n";
        }
#endif
        if (mps[mps.length()-1].col_dim()[0].first != right_end)
            throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    }
    
private:
    std::vector<int> basis_index;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
};

#endif
