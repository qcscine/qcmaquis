/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MPS_INITIALIZER_HELPER_H
#define MPS_INITIALIZER_HELPER_H

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/utils/BaseParameters.h"

#include <boost/tuple/tuple.hpp>

/**
 * @brief Helper class for the MPS initialization.
 *
 * This class wraps all the methods that are useful when initializing an MPS
 * from a given set of input data.
 *
 * For now, we include only a method, [GenerateIndexFromString], that converts
 * an input of integer (provided with the `init_basis_state` input parameter)
 * into a vector of tuples (charge, int). The charge is the symmetry block
 * which is populated upon construction, while the int is the position *within*
 * the symmetry block that is populated.
 *
 * By default, the [GenerateIndexFromString] method is deactivated.
 *
 * @tparam SymmGroup Symmetry group (for now we implement None and NU1 symmetry group)
 */
template<class SymmGroup>
class HelperClassBasisVectorConverter {
public:
  using ChargeType = typename SymmGroup::charge;
  using indexType = Index<SymmGroup>;
  using state_type = std::vector<boost::tuple<ChargeType, int> >;
  static state_type GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size)
  {
    throw std::runtime_error("GenerateIndexFromString method not available for this symmetry group");
  }
};

/** @brief Overload for the None class (to be used for vDMRG) */
template<>
class HelperClassBasisVectorConverter<TrivialGroup> {
public:
  // Types definition
  using indexType = Index<TrivialGroup>;
  using state_type = std::vector<boost::tuple<typename TrivialGroup::charge, int> >;
  // General implementation
  static state_type GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    auto state = state_type(size);
    for (int j = 0 ; j < size; ++j)
      state[j] = physDim[siteType[j]].element(inputVec[j]);
    return state;
  }
};

/** @brief Overload for the U1 class (to be used for vibronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<U1> {
public:
  // Types definition
  using indexType = Index<U1>;
  using state_type = std::vector<boost::tuple<typename U1::charge, int> >;
  // General implementation
  static state_type GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    auto state = state_type(size);
    // Note that here we have two possibilities: either a site is electronic site, or it is a vibrational one.
    for (int j = 0 ; j < size; ++j) {
      if (siteType[j] == 1) {
        auto posOfCharge = physDim[siteType[j]].position(inputVec[j]);
        state[j] = physDim[siteType[j]].element(posOfCharge);
      }
      else {
        state[j] = physDim[siteType[j]].element(inputVec[j]);
      }
    }
    return state;
  }
};

/**
 * @brief Overload of the previous class for the NU1 symmetry group.
 *
 * Note that, unlike in the previous case, where the input is given as a vector
 * of size L - L being the lattice size - here we give a vector of size
 * N - N being the template parameter for the NU1 class - and each element
 * is the position in the sublattice which is populated.
 *
 * @tparam N integer dimension of the NU1 class.
 */
template<int N>
class HelperClassBasisVectorConverter<NU1_template<N>> {
public:
  // Types definition
  using NU1 = NU1_template<N>;
  using indexType = Index<NU1>;
  using ChargeType = typename NU1::charge;
  using state_type = std::vector<boost::tuple<ChargeType, int> >;

  /** @brief Parser for the NU1 symmetry group
  * This function has two-fold functionality:
  * If the init_type is basis_state_generic, it returns the elements,
  * but if the init_type is basis_state_generic_const or basis_state_generic_default,
  * then it returns the charge and the integer to inidcate wether this site should be populated (1) or not (0)
  */
  static state_type GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != physDim.size())
      throw std::runtime_error("Index list number of elements does not match the number of site types. Check the setting 'init_basis_state'.");
    auto state = state_type(size);
    int numberOfTypes = inputVec.size();
    std::vector<int> counterOfTypes(numberOfTypes, 0);
    // If available, extracts the user-defined modals order
    std::vector<int> modalsOrder(size), inverseModalsOrder(size);
    if (!params.is_set("modals_order"))
      for (int p = 0; p < size; ++p)
        modalsOrder[p] = p;
    else
      modalsOrder = params["modals_order"].template as<std::vector<int> >();
    // Generates the inverse order
    for (int p = 0; p < modalsOrder.size(); ++p)
      inverseModalsOrder[p] = std::distance(modalsOrder.begin(), std::find(modalsOrder.begin(), modalsOrder.end(), p));
    // Fills the MPS.
    for (int iLattice = 0; iLattice < size; iLattice++) {
      auto positionOfSiteInNewLattice = inverseModalsOrder[iLattice];
      auto type = siteType[positionOfSiteInNewLattice];
      if (params["init_type"] == "basis_state_generic_const" || params["init_type"] == "basis_state_generic_default") {
        boost::tuple<ChargeType, bool> truePair = boost::make_tuple(boost::get<0>(physDim[type].element(0)), 1);
        boost::tuple<ChargeType, bool> falsePair = boost::make_tuple(boost::get<0>(physDim[type].element(0)), 0);   
        state[positionOfSiteInNewLattice] = (counterOfTypes[type] <= inputVec[type]) ? truePair : falsePair;
      } else {
        state[positionOfSiteInNewLattice] = (counterOfTypes[type] == inputVec[type]) ? physDim[type].element(0) : physDim[type].element(1);
      }
      counterOfTypes[type]++;
    }
    return state;
  }
};

#endif
