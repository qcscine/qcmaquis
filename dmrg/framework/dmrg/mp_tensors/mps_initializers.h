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
  // Declaration of types
  using StateType = std::vector<boost::tuple<typename SymmGroup::charge, int> >;

  /** @brief Class constructor */
  coherent_mps_init(BaseParameters & params, std::vector<Index<SymmGroup> > const& phys_dims_,
                    typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
    : coeff(params["init_coeff"].as<std::vector<double> >()), phys_dims(phys_dims_),
      site_type(site_type_), right_end(right_end_)
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
      state.resize(mps.length());
      for (int j=0; j<mps.length(); j++)
        state[j] = phys_dims[site_type[j]].element(basis_index[i][j]);
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
  }

private:
  typename SymmGroup::charge right_end;
  std::vector<double> coeff;
  std::vector<Index<SymmGroup> > phys_dims;
  std::vector<int> site_type;
  std::vector< std::vector<int> > basis_index;
  StateType state;
};

template<class Matrix, class SymmGroup>
class basis_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
  basis_mps_init(BaseParameters & params, std::vector<Index<SymmGroup> > const& phys_dims_,
                 std::vector<int> const& site_type_)
    : phys_dims(phys_dims_), site_type(site_type_)
  {
    std::string states = params["init_basis_state"].as<std::string>();
    std::vector<std::string> specifiedStates;
    boost::split(specifiedStates, states, boost::is_any_of("|"));
    std::stringstream ss(specifiedStates[0]);
    int ichar;
    while (ss >> ichar) {
        occupation.push_back(ichar);
        ss.ignore(1);
    }
  }


  void operator()(MPS<Matrix, SymmGroup> & mps)
  {
    assert(occupation.size() == mps.length());
    if (phys_dims[0].size() != 1)
      throw std::runtime_error("basis_mps_init only for TrivialGroup.");
    typedef typename SymmGroup::charge charge;
    charge C = SymmGroup::IdentityCharge;

    std::vector<boost::tuple<charge, int> > state(mps.length());
    for (int i=0; i<mps.length(); ++i)
        state[i] = boost::make_tuple(C, occupation[i]);
    mps = state_mps<Matrix>(state, phys_dims, site_type);
  }

private:
  std::vector<int> occupation;
  std::vector<Index<SymmGroup> > phys_dims;
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
    basis_mps_init_generic(BaseParameters & params_, const std::vector<Index<SymmGroup> >& phys_dims_,
                           typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
        : phys_dims(phys_dims_),
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
        mps = state_mps<Matrix>(state, phys_dims, site_type, right_end);
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
    BaseParameters& params;
};

template<class Matrix, class SymmGroup>
class basis_mps_init_generic_const : public mps_initializer<Matrix, SymmGroup>
{
  using state_type = std::vector<boost::tuple<typename SymmGroup::charge, int> >;

public:
  // -- Constructors --
  basis_mps_init_generic_const(BaseParameters & params_, const std::vector<Index<SymmGroup> >& phys_dims_,
                               typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
      : init_bond_dimension(params_["init_bond_dimension"]),
        phys_dims(phys_dims_), right_end(right_end_), site_type(site_type_), params(params_)
  {
    if (params["init_space"].str().empty())
      throw std::runtime_error("Init_space needs to be provided to populate basis_state_generic_const. Abort.");
    basis_index = params["init_space"].as<std::vector<int> >();
  }

  // Operator called when initialization occurs
  void operator()(MPS<Matrix, SymmGroup> & mps)
  {
    assert(basis_index.size() == mps.length());
    auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index, phys_dims, site_type, mps.length());
    // Actual MPS initialization
    mps = state_mps_const<Matrix>(state, phys_dims, site_type, right_end, false, init_bond_dimension);
    if (mps[mps.length()-1].col_dim()[0].first != right_end)
      throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    for (int i = 0; i < mps.length(); i++) {
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
  }
private:
  // -- ATTRIBUTES --
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
  // Types definition
  using state_type = std::vector<boost::tuple<typename SymmGroup::charge, int> >;

public:
  // -- Constructors --
  basis_mps_init_generic_default(BaseParameters & params_, std::vector<Index<SymmGroup> > const& phys_dims_,
                                 typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
      : init_bond_dimension(params_["init_bond_dimension"]),
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
    assert(basis_index.size() == mps.length());
    auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index, phys_dims, site_type, mps.length());
    mps = state_mps_const<Matrix>(state, phys_dims, site_type, right_end, true, init_bond_dimension);
    // Actual MPS initialization
    if (mps[mps.length()-1].col_dim()[0].first != right_end)
      throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    for (int i = 0; i < mps.length(); i++) {
      mps[i].divide_by_scalar(mps[i].scalar_norm());
    }
  }

private:
  // -- ATTRIBUTES --
  std::vector<int> basis_index;
  std::size_t init_bond_dimension;
  std::vector<Index<SymmGroup> > phys_dims;
  typename SymmGroup::charge right_end;
  std::vector<int> site_type;
  BaseParameters& params;
};

#endif
