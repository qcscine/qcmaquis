/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef MPS_INITIALIZER_POV_H
#define MPS_INITIALIZER_POV_H

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

template<class Matrix, class SymmGroup>
struct default_mps_init_pov : public mps_initializer_pov<Matrix, SymmGroup>
{
    // Types definition
    typedef typename MPS<Matrix, SymmGroup>::MPS MPSType ;
    typedef typename std::vector< MPSType >      MPSVector ;
    // Constructor
    default_mps_init_pov(BaseParameters & parms,
                         std::vector<Index<SymmGroup> > const& phys_dims_,
                         typename SymmGroup::charge right_end_,
                         std::vector<int> const& site_type_)
            : init_bond_dimension(parms["init_bond_dimension"])
            , phys_dims(phys_dims_)
            , right_end(right_end_)
            , site_type(site_type_)
    { }
    // Parenthesis operator
    void operator()(MPSVector & mps_vec)
    {
        init_sectors(mps_vec, this->init_bond_dimension, true);
    }
    // "Real" initialization routine
    void init_sectors(MPSVector & mps_vector,
                      size_t Mmax,
                      bool fillrand = true,
                      typename Matrix::value_type val = 0)
    {
        for (size_t idx=0; idx<mps_vector.size(); idx++) {
            std::size_t L = mps_vector[idx].length();
            parallel::scheduler_balanced scheduler(L);
            maquis::cout << "Right end: " << right_end << std::endl;
            std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, Mmax);
            omp_for(size_t i, parallel::range<size_t>(0, L), {
                parallel::guard proc(scheduler(i));
                mps_vector[idx][i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], fillrand, val);
                mps_vector[idx][i].divide_by_scalar(mps_vector[idx][i].scalar_norm());
            });
        }
    }
    // Attributes
    std::size_t init_bond_dimension;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
};

template<class Matrix, class SymmGroup>
class basis_mps_init_generic_pov : public mps_initializer_pov<Matrix, SymmGroup>
{
public:
    typedef std::vector<boost::tuple<typename SymmGroup::charge, size_t> > state_type;
    typedef typename MPS<Matrix, SymmGroup>::MPS MPSType ;
    typedef typename std::vector< MPSType >      MPSVector ;
    basis_mps_init_generic_pov(BaseParameters & params,
                              std::vector<Index<SymmGroup> > const& phys_dims_,
                              typename SymmGroup::charge right_end_,
                              std::vector<int> const& site_type_)
            : phys_dims(phys_dims_)
            , right_end(right_end_)
            , site_type(site_type_)
    {
        std::vector<std::string> list_sa;
        std::string input_str = params["follow_basis_state"].str();
        boost::split(list_sa, input_str, boost::is_any_of("|"));
        for (int i = 0; i < list_sa.size(); i++) {
            std::stringstream ss(list_sa[i]);
            int ichar;
            std::vector<int> tmp_vec;
            while (ss >> ichar) {
                tmp_vec.push_back(ichar);
                ss.ignore(1);
            }
            basis_index.push_back(tmp_vec);
        }
    }
    //
    // Operator used to actually initialize the vector
    // -----------------------------------------------
    void operator()(MPSVector & mps_vector)
    {
        for (int i = 0 ; i < basis_index.size() ; i++ ) {
            assert(basis_index[i].size() == mps_vector[i].length());
            state.resize(mps_vector[i].length());
            std::cout << "State average - state number " << i+1 << std::endl ;
            for (int j = 0 ; j < mps_vector[i].length() ; ++j) {
                state[j] = phys_dims[site_type[j]].element(basis_index[i][j]);
                maquis::cout << boost::get<0>(state[j]) << ":" << boost::get<1>(state[j])<< " ";
            }
            maquis::cout << "\n";
            mps_vector[i] = state_mps<Matrix>(state, phys_dims, site_type, right_end);
            if (mps_vector[i][mps_vector[i].length()-1].col_dim()[0].first != right_end)
                throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
        }
    }

private:
    std::vector< std::vector<int> > basis_index;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
    state_type state;
};

// +-------------------------+
//  INITIALIZER FROM CHK FILE
// +-------------------------+

template<class Matrix, class SymmGroup>
class chk_mps_init_pov : public mps_initializer_pov<Matrix, SymmGroup>
{
public:
    typedef typename MPS<Matrix, SymmGroup>::MPS MPSType ;
    typedef typename std::vector< MPSType >      MPSVector ;
    chk_mps_init_pov(BaseParameters & params)
    {
        boost::split(chk_names, params["follow_basis_state"].str(), boost::is_any_of("|"));
    }
    //
    // Operator used to actually initialize the vector
    // -----------------------------------------------
    void operator()(MPSVector & mps_vector)
    {
        assert (mps_vector.size() == chk_names.size()) ;
        for (std::size_t idx = 0; idx < mps_vector.size(); idx++)
            load(chk_names[idx], mps_vector[idx]) ;
    }

private:
    std::vector< std::string > chk_names ;
};


#endif
