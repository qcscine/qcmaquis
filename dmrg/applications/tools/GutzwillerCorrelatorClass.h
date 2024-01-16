/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Laboratory for Physical Chemistry, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@phys.chem.ethz.ch>
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

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

template <class Matrix, class SymmGroup>
class GutzwillerCalculator {
 public:
  // Types declaration
  using MPSType = MPS<Matrix, SymmGroup>;

  /**
   * @brief Class constructor
   * @param name1 hdf5 file where the first MPS is stored.
   * @param jValue correlation factor.
   */
  GutzwillerCalculator(std::string name, double jValue) : jValue_(jValue) {
    load(name, mps);
  }

  /** @brief Method to calculate the overlap */
  void applyCorrelator() {
    for (int iSite = 0; iSite < mps.size(); iSite++)
      mps[iSite].scaleByExponentialProductOfCharges(jValue_);
  }

  /** @brief Normalizes the underlying MPS wave function */
  void normalize() { mps[0] /= std::sqrt(norm(mps)); }

  /** @brief Getter for the MPS */
  auto getMPS() const { return mps; };

 private:
  MPSType mps;
  double jValue_;
};