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
  using state_type = std::vector<boost::tuple<ChargeType, std::size_t> >;
  static state_type GenerateIndexFromString(const std::vector<int>& inputVec, const std::vector<indexType>& physDim, 
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
  using state_type = std::vector<boost::tuple<typename TrivialGroup::charge, std::size_t> >;
  // General implementation
  static state_type GenerateIndexFromString(const std::vector<int>& inputVec, const std::vector<indexType>& physDim, 
                                            const std::vector<int>& siteType, int size) {
    assert(inputVec.size() == size);
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
  using state_type = std::vector<boost::tuple<typename U1::charge, std::size_t> >;
  // General implementation
  static state_type GenerateIndexFromString(const std::vector<int>& inputVec, const std::vector<indexType>& physDim, 
                                            const std::vector<int>& siteType, int size) {
    assert(inputVec.size() == size);
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
  using state_type = std::vector<boost::tuple<ChargeType, std::size_t> >;

  /** @brief Parser for the NU1 symmetry group */
  static state_type GenerateIndexFromString(const std::vector<int>& inputVec, const std::vector<indexType>& physDim, 
                                            const std::vector<int>& siteType, int size) {
    auto state = state_type(size);
    int numberOfTypes = inputVec.size();
    std::vector<int> counterOfTypes(numberOfTypes, 0);
    for (int iLattice = 0; iLattice < size; iLattice++) {
      if (counterOfTypes[siteType[iLattice]] == inputVec[siteType[iLattice]]) {
        state[iLattice] = physDim[siteType[iLattice]].element(0);
      }
      else {
        state[iLattice] = physDim[siteType[iLattice]].element(1);
      }
      counterOfTypes[siteType[iLattice]]++;
    }
    return state;
  }
};

#endif