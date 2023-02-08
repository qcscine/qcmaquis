/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiard <abaiard@ethz.ch>
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
