/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef MPS_INITIALIZER_HELPER_H
#define MPS_INITIALIZER_HELPER_H

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/utils/BaseParameters.h"

#include <boost/tuple/tuple.hpp>

namespace InitializerHelperFunctions {

/** @brief Method to convert a string into an ONV for the NU1 case */
template<int N>
auto GenerateIndexFromStringNMode(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<Index<NU1_template<N>>>& physDim,
                                  const std::vector<int>& siteType, int size)
{
  // Types definition
  using NU1 = NU1_template<N>;
  using indexType = Index<NU1>;
  using ChargeType = typename NU1::charge;
  using stateEntryType = std::vector<boost::tuple<ChargeType, int> >;
  using stateType = std::vector<stateEntryType>;
  // Main function body
  auto state = stateType(size, stateEntryType(1));
  if (inputVec.size() != physDim.size())
    throw std::runtime_error("Index list number of elements does not match the number of site types. Check the setting 'init_basis_state'.");
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
      state[positionOfSiteInNewLattice][0] = (counterOfTypes[type] <= inputVec[type]) ? truePair : falsePair;
    }
    else {
      state[positionOfSiteInNewLattice][0] = (counterOfTypes[type] == inputVec[type]) ? physDim[type].element(0) : physDim[type].element(1);
    }
    counterOfTypes[type]++;
  }
  return state;
}

} // namespace InitializerHelperFunctions

/**
 * @brief Helper class for the MPS initialization.
 *
 * This class wraps all the methods that are useful when initializing an MPS
 * from a given set of input data.
 *
 * For now, we include only a method, [GenerateIndexFromString], that converts
 * an input of integer (provided with the `init_basis_state` input parameter)
 * into a vector of vectors of tuples (charge, int). The charge is the symmetry block
 * which is populated upon construction, while the int is the position *within*
 * the symmetry block that is populated.
 *
 * By default, the [GenerateIndexFromString] method is deactivated.
 *
 * @tparam SymmGroup Symmetry group
 */
template<class SymmGroup>
class HelperClassBasisVectorConverter {
public:
  using ChargeType = typename SymmGroup::charge;
  using indexType = Index<SymmGroup>;
  using stateType = std::vector<std::vector<boost::tuple<ChargeType, int> > >;
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
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
  using stateEntryType = std::vector<boost::tuple<typename TrivialGroup::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    auto state = stateType(size, stateEntryType(1));
    std::cout << physDim[0] << std::endl;
    for (int j = 0 ; j < size; ++j) {
      state[j][0] = physDim[siteType[j]].element(inputVec[j]);
      std::cout << boost::get<0>(state[j][0]) << " " << boost::get<1>(state[j][0]) << std::endl;
    }
    return state;
  }
};

/** @brief Overload for the U1 class (to be used for vibronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<U1> {
public:
  // Types definition
  using indexType = Index<U1>;
  using stateEntryType = std::vector<boost::tuple<typename U1::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    auto state = stateType(size, stateEntryType(1));
    // Note that here we have two possibilities: either a site is electronic site, or it is a vibrational one.
    for (int j = 0 ; j < size; ++j) {
      if (siteType[j] == 1) {
        auto posOfCharge = physDim[siteType[j]].position(inputVec[j]);
        state[j][0] = physDim[siteType[j]].element(posOfCharge);
      }
      else {
        state[j][0] = physDim[siteType[j]].element(inputVec[j]);
      }
    }
    return state;
  }
};

/** @brief Overload for the TwoU1 class (to be used for electronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<TwoU1PG> {
public:
  // Types definition
  using indexType = Index<TwoU1PG>;
  using stateEntryType = std::vector<boost::tuple<typename TwoU1PG::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    std::vector<int> orbitalOrder(size);
    // Retrieves the orbital order
    if (!params.is_set("orbital_order"))
      for (int p = 0; p < size; ++p)
        orbitalOrder[p] = p+1;
    else
        orbitalOrder = params["orbital_order"].template as<std::vector<int> >();
    std::transform(orbitalOrder.begin(), orbitalOrder.end(), orbitalOrder.begin(), boost::lambda::_1-1);
    auto state = stateType(size, stateEntryType(1));
    for (int j = 0 ; j < size; ++j) {
      int hfIdx = inputVec[orbitalOrder[j]];
      if (hfIdx>4) {
        if (hfIdx==5 && (params["init_type"] == "basis_state_generic_const" || params["init_type"] == "basis_state_generic_default")) {
          state[j].resize(4);
          for (int addPhysDim = 0; addPhysDim < 4; ++addPhysDim)
            state[j][addPhysDim] = physDim[siteType[j]].element(addPhysDim);
        } else {
          throw std::runtime_error("HF coefficients range from 1 (empty)  to 4 (doubly occupied). 5 (mix of all possible occupations) only allowed for init_type generic_const/default");
        }
      } else {
        state[j][0] = physDim[siteType[j]].element(4-hfIdx);
      }
    }
    return state;
  }
};

/** @brief Overload for the SU2U1 class (to be used for electronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<SU2U1PG> {
public:
  // Types definition
  using indexType = Index<SU2U1PG>;
  using stateEntryType = std::vector<boost::tuple<typename SU2U1PG::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    std::vector<int> orbitalOrder(size);
    // Retrieves the orbital order
    if (!params.is_set("orbital_order"))
      for (int p = 0; p < size; ++p)
        orbitalOrder[p] = p+1;
    else
        orbitalOrder = params["orbital_order"].template as<std::vector<int> >();
    std::transform(orbitalOrder.begin(), orbitalOrder.end(), orbitalOrder.begin(), boost::lambda::_1-1);
    auto state = stateType(size, stateEntryType(1));
    for (int j = 0 ; j < size; ++j) {
      int hfIdx = inputVec[orbitalOrder[j]];
      maquis::cout << "In GenerateIndexFromString with hfIdx: " << hfIdx << std::endl;
      if (hfIdx>4) {
        if (hfIdx==5 && (params["init_type"] == "basis_state_generic_const" || params["init_type"] == "basis_state_generic_default")) {
          state[j].resize(4);
          for (int addPhysDim = 0; addPhysDim < 4; ++addPhysDim)
            state[j][addPhysDim] = physDim[siteType[j]].element(addPhysDim);
        } else if (hfIdx==6 || hfIdx==7) { // to explicitly construct a specified csf
          state[j].resize(1);
          state[j][0] = (hfIdx==6) ? physDim[siteType[j]].element(1) : physDim[siteType[j]].element(2);
        } else {
          throw std::runtime_error("HF coefficients range from 1 (empty)  to 4 (doubly occupied). 5 (mix of all possible occupations) only allowed for init_type generic_const/default");
        }
      } else if (hfIdx==2 || hfIdx==3) { //since in this case alpha and beta are equivalent
        state[j].resize(2);
        state[j][0] = physDim[siteType[j]].element(1);
        state[j][1] = physDim[siteType[j]].element(2);
      } else { //otherwise the occupation is either zero or double
        state[j].resize(1);
        state[j][0] = physDim[siteType[j]].element(4-hfIdx);
      }
    }
    return state;
  }
};

/** @brief Overload for the SU2U1 class (to be used for electronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<SU2U1> {
public:
  // Types definition
  using indexType = Index<SU2U1>;
  using stateEntryType = std::vector<boost::tuple<typename SU2U1::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size) {
    if (inputVec.size() != size)
      throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
    std::vector<int> orbitalOrder(size);
    // Retrieves the orbital order
    if (!params.is_set("orbital_order"))
      for (int p = 0; p < size; ++p)
        orbitalOrder[p] = p+1;
    else
        orbitalOrder = params["orbital_order"].template as<std::vector<int> >();
    std::transform(orbitalOrder.begin(), orbitalOrder.end(), orbitalOrder.begin(), boost::lambda::_1-1);
    auto state = stateType(size, stateEntryType(1));
    for (int j = 0 ; j < size; ++j) {
      int hfIdx = inputVec[orbitalOrder[j]];
      if (hfIdx>4) {
        if (hfIdx==5 && (params["init_type"] == "basis_state_generic_const" || params["init_type"] == "basis_state_generic_default")) {
          state[j].resize(4);
          for (int addPhysDim = 0; addPhysDim < 4; ++addPhysDim)
            state[j][addPhysDim] = physDim[siteType[j]].element(addPhysDim);
        } else if (hfIdx==6 || hfIdx==7) { // to explicitly construct a specified csf
          state[j].resize(1);
          state[j][0] = (hfIdx==6) ? physDim[siteType[j]].element(1) : physDim[siteType[j]].element(2);
        } else {
          throw std::runtime_error("HF coefficients range from 1 (empty)  to 4 (doubly occupied). 5 (mix of all possible occupations) only allowed for init_type generic_const/default");
        }
      } else if (hfIdx==2 || hfIdx==3) { //since in this case alpha and beta are equivalent
        state[j].resize(2);
        state[j][0] = physDim[siteType[j]].element(1);
        state[j][1] = physDim[siteType[j]].element(2);
      } else {
        state[j].resize(1);
        state[j][0] = physDim[siteType[j]].element(4-hfIdx);
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
  using stateEntryType = std::vector<boost::tuple<ChargeType, int> >;
  using stateType = std::vector<stateEntryType>;

  /** @brief Parser for the NU1 symmetry group
  * This function has two-fold functionality:
  * If the init_type is basis_state_generic, it returns the elements,
  * but if the init_type is basis_state_generic_const or basis_state_generic_default,
  * then it returns the charge and the integer to inidcate wether this site should be populated (1) or not (0)
  */
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size)
  {
    return InitializerHelperFunctions::GenerateIndexFromStringNMode<N>(params, inputVec, physDim, siteType, size);
  }
};

/** @brief Overload for the TwoU1 class (to be used for electronic Hamiltonians) */
template<>
class HelperClassBasisVectorConverter<TwoU1> {
public:
  // Types definition
  using indexType = Index<TwoU1>;
  using stateEntryType = std::vector<boost::tuple<typename TwoU1::charge, int> >;
  using stateType = std::vector<stateEntryType>;
  using ChargeType = typename TwoU1::charge;
  // General implementation
  static stateType GenerateIndexFromString(BaseParameters& params, const std::vector<int>& inputVec, const std::vector<indexType>& physDim,
                                            const std::vector<int>& siteType, int size)
  {
    // Return object definition
    auto state = stateType(size, stateEntryType(1));
    // We have to be careful in the TwoU1 case, because we the case in which TwoU1 is the electronic Hamiltonian, we have
    // to invoke a specialized constructor. For PreBO and nMode vDMRG, we can instead invoke the general function
    if (params["MODEL"] == "quantum_chemistry") {
      if (inputVec.size() != size)
        throw std::runtime_error("Index list number of elements does not match the lattice size. Check the input settings.");
      std::vector<int> orbitalOrder(size);
      // Retrieves the orbital order
      if (!params.is_set("orbital_order"))
        for (int p = 0; p < size; ++p)
          orbitalOrder[p] = p+1;
      else
          orbitalOrder = params["orbital_order"].template as<std::vector<int> >();
      std::transform(orbitalOrder.begin(), orbitalOrder.end(), orbitalOrder.begin(), boost::lambda::_1-1);
      for (int j = 0 ; j < size; ++j) {
        int hfIdx = inputVec[orbitalOrder[j]];
        if (hfIdx > 4) {
          if (hfIdx==5 && (params["init_type"] == "basis_state_generic_const" || params["init_type"] == "basis_state_generic_default")) {
            state[j].resize(4);
            for (int addPhysDim = 0; addPhysDim < 4; ++addPhysDim)
              state[j][addPhysDim] = physDim[siteType[j]].element(addPhysDim);
          }
          else {
            throw std::runtime_error("HF coefficients range from 1 (empty)  to 4 (doubly occupied). 5 (mix of all possible occupations) only allowed for init_type generic_const/default");
          }
        }
        else {
          state[j][0] = physDim[siteType[j]].element(4-hfIdx);
 