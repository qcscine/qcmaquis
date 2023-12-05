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
  public:
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
      std::size_t L = mps.length();
      for (int i = 0; i < L; i++) {
        parallel::guard proc(scheduler(i));
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], fillrand, val);
        mps[i].divide_by_scalar(mps[i].scalar_norm());
      }
    }

    // Class members
    int init_bond_dimension;
    std::vector<Index<SymmGroup> > phys_dims;
    typename SymmGroup::charge right_end;
    std::vector<int> site_type;
};

/**
 * @brief Const MPS initializer (MPS with all entries = 1.) 
 */
template<class Matrix, class SymmGroup>
struct const_mps_init : public mps_initializer<Matrix, SymmGroup>
{
  public:
    /** @brief Class costructor */
    const_mps_init(BaseParameters & parms, std::vector<Index<SymmGroup> > const& phys_dims,
                  typename SymmGroup::charge right_end, std::vector<int> const& site_type)
      : di(parms, phys_dims, right_end, site_type) { }

    void operator()(MPS<Matrix, SymmGroup> & mps) { di.init_sectors(mps, di.init_bond_dimension, false, 1.); }

    // Class member
    default_mps_init<Matrix, SymmGroup> di;
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
      std::string states = params["init_basis_state"].template as<std::string>();
      std::vector<std::string> specifiedStates;
      boost::split(specifiedStates, states, boost::is_any_of("|"));
      std::stringstream ss(specifiedStates[0]);
      int ichar;
      if (params["init_state_type"] == "csf" && (sym=="su2u1" || sym=="su2u1pg")){
        char jchar;
        while (ss.get(jchar)) {
          if (jchar == '2') ichar = 4; // doubly occ to explicitly construct csf
          else if (jchar == 'u') ichar = 6; // up to explicitly construct csf
          else if (jchar == 'd') ichar = 7; // down to explicitly construct csf
          else if  (jchar == '0') ichar = 1; // not occ to explicitly construct csf
          else throw std::runtime_error("The specificed state contains symbols which are not recognized. Abort.");
          basis_index.push_back(ichar);
          ss.ignore(1);
        }
      } else { // regular determinant
        while (ss >> ichar) {
          basis_index.push_back(ichar);
          ss.ignore(1);
        }
      }
    }


    /** @brief Operator (), called when the MPS is constructed */
    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
      assert(basis_index.size() == mps.length());
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
    basis_index = params["init_space"].template as<std::vector<int> >();
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
    basis_index = params["init_space"].template as<std::vector<int> >();
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


/** @brief Coherent MPS initialization */
template<class Matrix, class SymmGroup>
class coherent_mps_init : public mps_initializer<Matrix, SymmGroup>
{
  public:
    /** @brief Class constructor 
     *
     * Note that generally determinants are utilized except for SU2, where CSF are constructed.
     * 
     */
    coherent_mps_init(BaseParameters & params_, std::vector<Index<SymmGroup> > const& phys_dims_,
                      typename SymmGroup::charge right_end_, std::vector<int> const& site_type_)
      : phys_dims(phys_dims_), site_type(site_type_), right_end(right_end_), params(params_)
    {
      if (params["init_file"].str().empty() && params["init_basis_state"].str().empty())
        throw std::runtime_error("Either init_file or init_basis_state has to be provided for the coherent initializer");

      initialBondDim = (params["init_bond_dimension"] > 5) ? params["init_bond_dimension"] : params["max_bond_dimension"];
      fromFile = (params["init_file"].str().empty()) ? false : true;

      std::vector<std::string> list_dets;

      if (fromFile) {
        std::string fileName = params["init_file"];
        std::vector<std::string> specifiedFiles;
        boost::split(specifiedFiles, fileName, boost::is_any_of("|"));
        fileName=specifiedFiles[0]; // This is a safeguard for the interface sim initialization in case several filenames are provided
        if (!boost::filesystem::exists(fileName))
          throw std::runtime_error("Initializer file " + fileName + " does not exist\n");
        maquis::cout << "Initializing MPS from file " << fileName << std::endl;
        std::ifstream stateFile;
        stateFile.open(fileName.c_str());
        std::string line;
        std::vector< std::string > line_splitted;
        while (std::getline(stateFile, line)) {
          boost::trim_left(line);
          boost::trim_right(line);
          boost::split(line_splitted, line, boost::is_any_of(" "), boost::token_compress_on);
          coeffs.push_back(std::stod(line_splitted[0]));
          list_dets.push_back(line_splitted[1]);
        }
        stateFile.close();
      } else {
        boost::split(list_dets, params["init_basis_state"].str(), boost::is_any_of("|"));
        coeffs = params_["init_coeffs"].as<std::vector<double> >();
      }

      for (int i = 0; i < list_dets.size(); i++) {
        std::stringstream ss(list_dets[i]);
        int ichar;
        std::vector<int> tmp_vec;
        if (params["init_state_type"] == "csf"){
          char jchar;
          while (ss.get(jchar)) {
            if (jchar == '2') ichar = 4; // doubly occ to explicitly construct csf
            else if (jchar == 'u') ichar = 6; // up to explicitly construct csf
            else if (jchar == 'd') ichar = 7; // down to explicitly construct csf
            else if  (jchar == '0') ichar = 1; // not occ to explicitly construct csf
            else throw std::runtime_error("The specificed state contains symbols which are not recognized. Abort.");
            tmp_vec.push_back(ichar);
            ss.ignore(1);
          }
        } else { // regular determinant
          while (ss >> ichar) {
            tmp_vec.push_back(ichar);
            ss.ignore(1);
          }
        }
        basis_index.push_back(tmp_vec);
      }
      // Final check
      assert (basis_index.size() == coeffs.size());
    }

    /** @brief Method to construct the MPS */
    void operator()(MPS<Matrix, SymmGroup>& mps)
    {
      MPS<Matrix, SymmGroup> MPSBuffer;
      auto sym = params["symmetry"].as<std::string>();
      for (int i=0; i<basis_index.size(); i++ ) {
        auto state = HelperClassBasisVectorConverter<SymmGroup>::GenerateIndexFromString(params, basis_index[i], phys_dims, site_type, mps.length());
        auto mps_tmp = (sym=="su2u1" || sym=="su2u1pg") ? state_mps_cd<Matrix>(state, phys_dims, site_type, right_end, initialBondDim, false)
                                                        : state_mps<Matrix>(state, phys_dims, site_type, right_end);
        mps_tmp.normalize_right();
        if (i == 0) {
          mps = mps_tmp;
          mps[0] *= coeffs[0];
        }
        else {
          MPSBuffer = join(mps, mps_tmp, 1., coeffs[i]);
          mps = MPSBuffer;
        }
      }
      // Compression to initialBondDim
      mps = compression::l2r_compress(mps, initialBondDim, 0);
      // Normalization
      for (int i = 0; i < mps.length(); i++) {
        mps[i].divide_by_scalar(mps[i].scalar_norm());
      }
      if (mps[mps.length()-1].col_dim()[0].first != right_end)
        throw std::runtime_error("Initial state does not satisfy total quantum numbers.");
    }

  private:
    typename SymmGroup::charge right_end;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_type;
    std::vector< std::vector<int> > basis_index;
    std::vector<double> coeffs;
    BaseParameters& params;
    int initialBondDim;
    bool fromFile;
};

#endif
