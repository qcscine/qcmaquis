/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_NMODE_LATTICE
#define MAQUIS_DMRG_NMODE_LATTICE

#ifdef DMRG_VIBRATIONAL

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"
#include "LatticeHelperClass.hpp"

/**
 * @brief Lattice representing a n-mode vibrational Hamiltonian.
 *
 * Tha lattice has the following properties:
 *  - each site is mapped to a modal which is, in turn, associated with a mode.
 *  - the lattice is, therefore, partition in sublattices, one associated with
 *    each mode.
 *  - in the conventional sorting, modals associated with the same mode are
 *    close one to each other. However, the sorting can be modified such that
 *    modals associated with different modes are "intertwined".
 */

class NModeLattice : public lattice_impl
{
public:
  // Types definition
  using post_t = lattice_impl::pos_t;

  /**
   * @brief Class constructor for the lattice
   * @param parameters parameter container
   */
  NModeLattice(BaseParameters& parameters) : L(0), vectorWithStartingPositions(0), vector_types(0)
  {
    // Parameter initialization
    L = parameters["L"];
    order = LatticeHelperClass::getOrbitalOrder(parameters, "modals_order", false);
    inverseOrder.resize(L);
    for (int iOrder = 0; iOrder < order.size(); iOrder++) {
      auto pointerToFoundElement = std::find(order.begin(), order.end(), iOrder);
      assert(pointerToFoundElement != order.end());
      inverseOrder[iOrder] = std::distance(order.begin(), pointerToFoundElement);
    }
    num_modes = parameters["nmode_num_modes"];
    std::string jnk = parameters["nmode_num_basis"];
    numTypes = num_modes;
    // The vector [size_vec] contains the number of modal bases per mode.
    std::vector<std::string> size_vec;
    boost::split(size_vec, jnk, boost::is_any_of(","));
    assert(size_vec.size() == num_modes);
    // Loops over the number of basis functions and calculates the index
    // of the first site associated with a given mode
    vectorWithStartingPositions.reserve(L);
    for (int idx = 0; idx < size_vec.size(); idx++)
    {
      if (idx == 0) {
          vectorWithStartingPositions.push_back(0);
      }
      else {
        int mod = stoi(size_vec[idx-1]);
        if (mod <= 0)
          throw std::runtime_error("Non-positive number of basis function found");
        else
          vectorWithStartingPositions.push_back(vectorWithStartingPositions[idx-1] + mod);
      }
    }
    int count = vectorWithStartingPositions[size_vec.size()-1] + stoi(size_vec[size_vec.size()-1]);
    if (count != L)
      throw std::runtime_error("Inconsistent number of basis functions");
    // Now populates the vector with the type of each site (i.e., the mode to
    // which they belong
    vector_types.resize(L);
    int jcont=0;
    // Loop over the modes
    int iMode = 0;
    for (int iSite = 0; iSite < L; iSite++) {
        if (iMode != num_modes-1 && iSite == vectorWithStartingPositions[iMode+1])
            iMode++;
        vector_types[inverseOrder[iSite]] = iMode;
    }
  }

  /** @brief Returns the next position in the lattice */
  std::vector<pos_t> forward(pos_t i) const override {
    std::vector<pos_t> ret;
    if (i < L-1)
        ret.push_back(i+1);
    return ret;
  }

  /** @brief Returns the neighbors of a given site */
  std::vector<pos_t> all(pos_t i) const override {
    std::vector<pos_t> ret;
    if (i < L-1)
        ret.push_back(i+1);
    if (i > 0)
        ret.push_back(i-1);
    return ret;
  }

  /**
   * @brief Getter for the property
   *
   * Note that, in additional to the usual properties of a lattice, we code
   * the additional property "sublatticePos" which states where the sublattice
   * associated with a given mode is starting.
   * Also, we define a property called "absolutePositionInLattice" that, given
   * the mode and modal index, returns the corresponding absolute position in the
   * DMRG lattice.
   *
   * @param property string identifier for the property
   * @param pos vector of positions
   * @return boost::any requested property
   */
  boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const override
  {
    if (property == "label" && pos.size() == 1)
      return boost::any(site_label(pos[0]));
    else if (property == "label" && pos.size() == 2)
      return boost::any(bond_label(pos[0], pos[1]));
    else if (property == "type" && pos.size() == 1)
      return boost::any(vector_types[pos[0]]);
    else if (property == "type" && pos.size() == 2)
      return boost::any(0);
    else if (property == "absolutePositionInLattice" && pos.size() == 2) {
      int posInConventionalSorting = vectorWithStartingPositions[pos[0]] + pos[1];
      return boost::any(inverseOrder[posInConventionalSorting]);
    }
    else if (property == "sublatticePos" && pos.size() == 1) {
      assert (pos[0] >= 0 && pos[0] < L);
      return boost::any(vectorWithStartingPositions[pos[0]]);
    }
    else if (property == "ParticleType" && pos.size() == 1) {
      assert (pos[0] >= 0 && pos[0] < L);
      return vector_types[pos[0]];
    }
    else if (property == "NumTypes") {
      return num_modes;
    }
    else {
      std::ostringstream ss;
      ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
      throw std::runtime_error(ss.str());
      return boost::any();
    }
  }

  /** @brief Getter for the lattice size */
  pos_t size() const override { return L; }

  /** @brief Getter for the number of types of sites */
  int getMaxType() const override { return numTypes; }

private:
  /** Size of the DMRG lattice */
  pos_t L;
  /** Largest index for site types */
  int numTypes;
  /** Number of modes (== number of site types) */
  int num_modes;
  /** Order of the modals in the DMRG lattice */
  std::vector<int> order, inverseOrder;
  /** Sites type vector */
  std::vector<int> vector_types;
  /**
   * The i-th element returns the number of basis that have been used
   * before the i-th mode. Used as offset in vectors
   */
  std::vector<int> vectorWithStartingPositions;

  /** @brief Prints the label of a given site */
  std::string site_label (int i) const {
      return "( " + boost::lexical_cast<std::string>(i) + " )";
  }

  /** @brief Prints the label of a given bond */
  std::string bond_label (int i, int j) const {
      return (  "( " + boost::lexical_cast<std::string>(i) + " )"
              + " -- "
              + "( " + boost::lexical_cast<std::string>(j) + " )");
  }
};

#endif // DMRG_VIBRATIONAL

#endif // MAQUIS_DMRG_NMODE_LATTICE
