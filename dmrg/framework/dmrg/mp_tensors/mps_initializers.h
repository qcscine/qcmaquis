/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
  coherent_mps_init(BaseParameters & params, std::vector<Index<SymmGroup> > const& phys_dims_,
                    std::vector<int> const& site_type_)
    : coeff(params["init_coeff"].as<std::vector<double> >()), phys_dims(phys_dims_),
      site_type(site_type_) { }

  void operator()(MPS<Matrix, SymmGroup> & mps)
  {
    // Types definition
    typedef typename SymmGroup::charge charge;
    using std::exp;
    using std::sqrt;
    using std::pow;
    using boost::math::factorial;
    // Checks
    assert(coeff.size() == mps.length());
    if (phys_dims[0].size() != 1)
      throw std::runtime_error("coherent_mps_init only for TrivialGroup.");
    // Variable extraction
    auto L = coeff.size();
    Index<SymmGroup> trivial_i;
    trivial_i.insert(std::make_pair(SymmGroup::IdentityCharge, 1));
    // MPS Initialization
    for (int p=0; p<L; ++p) {
      int s=0;
      Matrix m(phys_dims[site_type[p]][s].second, 1, 0.);
      for (int ss=0; ss<phys_dims[site_type[p]][s].second; ++ss)
        m(ss, 0) = pow(coeff[p], ss) * sqrt(factorial<double>(ss)) / factorial<double>(ss);
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

/**
 * @brief ONV MPS initializer
 *
 * The MPS is initialized from a unique ONV.
 * The initial MPS has, therefore, m=1
 */
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
