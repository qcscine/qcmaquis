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
#include "dmrg/utils/random.hpp"

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
    : init_bond_dimension(parms["init_bond_dimension"]), phys_dims(phys_dims_), right_end(right_end_),
      site_type(site_type_)
  {
    if (parms.is_set("seed"))
      dmrg_random::engine.seed(parms["seed"]);
  }

  /**
   * @brief Functor operator called to generate the MPS
   *
   * Note that fillrand variable is set here to true - so random initialization is done
   * by default.
   *
   * @param mps output MPS
   */
  void operator()(MPS<Matrix, SymmGroup> & mps) { init_sectors(mps, this->init_bond_dimension, true); }

  // Main routine
  void init_sectors(MPS<Matrix, SymmGroup> & mps, size_t Mmax, bool fillrand=true, typename Matrix::value_type val=0)
  {
    parallel::scheduler_balanced scheduler(mps.length());
    // Compute the indexes which are allowed by symmetry
    std::vector<Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, Mmax);
    // Populates the MPS tensor
    //omp_for(size_t i, parallel::range<size_t>(0, L), {
    std::size_t L = mps.length();
    for (int i = 0; i < L; i++) {
      parallel::guard proc(scheduler(i));
      mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], fillrand, val);
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
    //});
  }

  // Class members
  int init_bond_dimension;
  std::vector<Index<SymmGroup> > phys_dims;
  typename SymmGroup::charge right_end;
  std::vector<int> site_type;
};

/** @brief Const MPS initializer (MPS with all entries = 1.) */
template<class Matrix, class SymmGroup>
struct const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
  /** @brief Class costructor */
  const_mps_init(BaseParameters & parms, std::vector<Index<SymmGroup> > const& phys_dims,
                 typename SymmGroup::charge right_end, std::vector<int> const& site_type)
    : di(parms, phys_dims, right_end, site_type) { }

  void operator()(MPS<Matrix, SymmGroup> & mps) { di.init_sectors(mps, di.init_bond_dimension, false, 1.); }

  // Class member
  default_mps_init<Matrix, SymmGroup> di;
};

/** @brief Coherent MPS initialization */
template<class Matrix, class SymmGroup>
class coherent_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
  /** @brief Class constructor */
  coherent_mps_init(BaseParameters & params_, std::vector<Index<SymmGroup> > const& phys_dims_,
                    typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
    : coeff(params["init_coeff"].as<std::vector<double> >()), phys_dims(phys_dims_),
      site_type(site_type_), right_end(right_end_), params(params_)
  {
    std::vector<std::string> list_sa;
    std::string input_str = params["init_basis_state"].str();
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
    // Final check
    assert (basis_index.size() == coeff.size());
  }

  /** @brief Method to construct the MPS */
  void operator()(MPS<Matrix, SymmGroup>& mps)
  {
    MPS<Matrix, SymmGroup> MPSBuffer;
    for (int i=0; i<basis_index.size(); i++ ) {
      auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index[i], phys_dims, site_type, mps.length());
      auto mps_tmp = state_mps<Matrix>(state, phys_dims, site_type, right_end);
      if (i == 0) {
        mps = mps_tmp;
        mps[0] *= coeff[0];
      }
      else {
        MPSBuffer = join(mps, mps_tmp, 1., coeff[i]);
        mps = MPSBuffer;
      }
    }
    if (mps[mps.length()-1].col_dim()[0].first != right_end)
      throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    for (int i = 0; i < mps.length(); i++) {
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
  }

private:
  typename SymmGroup::charge right_end;
  std::vector<double> coeff;
  std::vector<Index<SymmGroup> > phys_dims;
  std::vector<int> site_type;
  std::vector< std::vector<int> > basis_index;
  BaseParameters& params;
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
    /**
     * @brief Class constructor from a parameter object
     * @param params Parameter container.
     * @param phys_dims_ Vector with the physical index per site type.
     * @param right_end_ Overall symmetry sector to which the MPS belongs.
     * @param site_type_ Vector with size == the lattice size, with the type of each site.
     */
    basis_mps_init_generic(BaseParameters & params_, const std::vector<Index<SymmGroup> >& phys_dims_,
                           typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
        : sym(params_["symmetry"].as<std::string>()), phys_dims(phys_dims_),
          right_end(right_end_), site_type(site_type_), params(params_)
    {
      std::string states = params["init_basis_state"].as<std::string>();
      std::vector<std::string> specifiedStates;
      boost::split(specifiedStates, states, boost::is_any_of("|"));
      std::stringstream ss(specifiedStates[0]);
      int ichar;
      while (ss >> ichar) {
        basis_index.push_back(ichar);
        ss.ignore(1);
      }
    }


    /** @brief Operator (), called when the MPS is constructed */
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        // assert(basis_index.size() == mps.length());
        auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index, phys_dims, site_type, mps.length());
        if (sym=="su2u1" || sym=="su2u1pg") { // SU2 electronic case --> special because of spin symmetries etc --> directly use state_mps_cd
          mps = state_mps_cd<Matrix>(state, phys_dims, site_type, right_end, 1, false);
        } else {
          mps = state_mps<Matrix>(state, phys_dims, site_type, right_end, 1);
        }
        if (mps[mps.length()-1].col_dim()[0].first != right_end)
            throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
        
        for (int i = 0; i < mps.length(); i++) {
            mps[i].divide_by_scalar(mps[i].scalar_norm());
        }
    }

private:
    std::string sym;
    std::vector<int> basis_index;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
    BaseParameters& params;
};

template<class Matrix, class SymmGroup>
class basis_mps_init_generic_const : public mps_initializer<Matrix, SymmGroup>
{
public:
  // -- Constructors --
  basis_mps_init_generic_const(BaseParameters & params_, const std::vector<Index<SymmGroup> >& phys_dims_,
                               typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
      : sym(params_["symmetry"].as<std::string>()), init_bond_dimension(params_["init_bond_dimension"]),
        phys_dims(phys_dims_), right_end(right_end_), site_type(site_type_), params(params_)
  {
    if (params["init_space"].str().empty())
      throw std::runtime_error("Init_space needs to be provided to populate basis_state_generic_const. Abort.");
    basis_index = params["init_space"].as<std::vector<int> >();
  }

  // Operator called when initialization occurs
  void operator()(MPS<Matrix, SymmGroup> & mps)
  {
    auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index, phys_dims, site_type, mps.length());
    assert(state.size() == mps.length());
    // Actual MPS initialization
    if (sym=="2u1" || sym=="2u1pg" || sym=="su2u1" || sym=="su2u1pg") { // electronic case --> special because of spin symmetries etc --> directly use state_mps
      mps = state_mps_cd<Matrix>(state, phys_dims, site_type, right_end, init_bond_dimension, false);
    } else {
      mps = state_mps_const<Matrix>(state, phys_dims, site_type, right_end, false, init_bond_dimension);
    }
    if (mps[mps.length()-1].col_dim()[0].first != right_end)
      throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    for (int i = 0; i < mps.length(); i++) {
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
  }
private:
  // -- ATTRIBUTES --
  std::string sym;
  std::vector<int> basis_index;
  std::size_t init_bond_dimension;
  std::vector<Index<SymmGroup> > phys_dims;
  typename SymmGroup::charge right_end;
  std::vector<int> site_type;
  BaseParameters& params;
};

template<class Matrix, class SymmGroup>
class basis_mps_init_generic_default : public mps_initializer<Matrix, SymmGroup>
{
public:
  // -- Constructors --
  basis_mps_init_generic_default(BaseParameters & params_, std::vector<Index<SymmGroup> > const& phys_dims_,
                                 typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
      : sym(params_["symmetry"].as<std::string>()), init_bond_dimension(params_["init_bond_dimension"]),
        phys_dims(phys_dims_), right_end(right_end_), site_type(site_type_), params(params_)
  {
    if (params["init_space"].str().empty())
      throw std::runtime_error("Init_space needs to be provided to populate basis_state_generic_default. Abort.");
    basis_index = params["init_space"].as<std::vector<int> >();
    if (params.is_set("seed"))
      dmrg_random::engine.seed(params["seed"]);
  }
  // Operator called when initialization occurs
  void operator()(MPS<Matrix, SymmGroup> & mps)
  {
    auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index, phys_dims, site_type, mps.length());
    assert(state.size() == mps.length());
    // Actual MPS initialization
    if (sym=="2u1" || sym=="2u1pg" || sym=="su2u1" || sym=="su2u1pg") { // electronic case --> special because of spin symmetries etc --> directly use state_mps
      mps = state_mps_cd<Matrix>(state, phys_dims, site_type, right_end, init_bond_dimension, true);
    } else {
      mps = state_mps_const<Matrix>(state, phys_dims, site_type, right_end, true, init_bond_dimension);
    }
    if (mps[mps.length()-1].col_dim()[0].first != right_end)
      throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    for (int i = 0; i < mps.length(); i++) {
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
  }

private:
  // -- ATTRIBUTES --
  std::string sym;
  std::vector<int> basis_index;
  std::size_t init_bond_dimension;
  std::vector<Index<SymmGroup> > phys_dims;
  typename SymmGroup::charge right_end;
  std::vector<int> site_type;
  BaseParameters& params;
};

#endif
