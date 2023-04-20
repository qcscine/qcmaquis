/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef LATTICE_HELPER_CLASS
#define LATTICE_HELPER_CLASS

#include "dmrg/utils/BaseParameters.h"

/** @brief Helper class with the functions used by all lattices */
class LatticeHelperClass {
public:

  /**
   * @brief Extracts the orbital order from a parameter container.
   * @param parms Parameter container object.
   * @param nameOfOrder Name of the parameter where to look for the order.
   * @param doShift if true, shift by -1 the order (this is needed if the order starts from 0)
   * @return std::vector<std::size_t> Order specificed in the input
   */
  static auto getOrbitalOrder(BaseParameters& parms, std::string nameOfOrder, bool doShift)
  {
    // Definition of key parameters
    using PositionType = int;
    int latticeSize = parms["L"];
    std::vector<PositionType> outputOrder(latticeSize, 0);
    // Main code part
    if (!parms.is_set(nameOfOrder)) {
      for (int p = 0; p < latticeSize; ++p)
        outputOrder[p] = p;
    }
    else {
      outputOrder = parms[nameOfOrder].as<std::vector<PositionType> >();
      if (outputOrder.size() != latticeSize)
        throw std::runtime_error("Number of orbitals in the orbital order does not match the total number of orbitals");
      // Shifts the order by -1 to match the convention that indices start in
      // C++ from *zero*.
      if (doShift)
        for (auto&& o: outputOrder)
          o--;
    }
    return outputOrder;
  }
};

#endif
