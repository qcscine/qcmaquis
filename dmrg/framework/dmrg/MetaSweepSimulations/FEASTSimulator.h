/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef FEAST_SIMULATOR
#define FEAST_SIMULATOR

#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/utils/BaseParameters.h"

template<class Matrix, class SymmGroup>
class FEASTSimulator {
  using LatticeType = Lattice;
  using ModelType = Model<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
  using MPOType = MPO<Matrix, SymmGroup>;
  using ValueType = typename MPSType::value_type;
  using InitializerType = mps_initializer<Matrix, SymmGroup>;

public:

  /** @brief Class constructor */
  FEASTSimulator(BaseParameters& parms, const ModelType& model, const LatticeType& lattice) : currentIter(0) {
    // Retrieve simulation parameters
    numStates = parms["feast_num_states"].as<int>();
    maxFeastIter = parms["feast_max_iter"].as<int>();
    eMin = parms["feast_emin"].as<double>();
    eMax = parms["feast_emax"].as<double>();
    numQuadraturePoint = parms["feast_num_points"].as<int>();
    intModality = parms["feast_integral_type"].as<std::string>();
    truncModality = parms["feast_truncation_type"].as<std::string>();
    // Checks consistency of the input
    if (intModality != "half" && intModality != "full")
      throw std::runtime_error("Parameter [feast_integral_type] not recognized");
    if (truncModality != "each" && truncModality != "end")
      throw std::runtime_error("Parameter [feast_truncation_type] not recognized");
    // Generates the initial guess for the MPSs
    initializeGuess(parms);
  }

private:

  /** @brief Generates the guess for FEAST */
  void initializeGuess(BaseParameters& parms) {
  }


  // Class members
  int currentIter;                // Index of the current FEAST iteration.
  int numStates;                  // Number of states to be targeted.
  int maxFeastIter;               // Maximum number of FEAST iterations.
  double eMin;                    // Lower bound for the complex contour integral
  double eMax;                    // Upper bound for the complex contour integral
  double feastThreshold;          // Threshold to assess the convergence of DMRG[FEAST]
  int numQuadraturePoint;         // Number of quadrature point
  std::string intModality;        // "Full" for the full circle integration, "half" for the half-circle one
  std::string truncModality;      // "Each" if the MPS must be truncated after each sum, "end" if the truncation must be done only at the end
  std::vector<MPSType> mpsGuess;  // Stores the current guess for hte FEAST procedure
};

#endif // FEAST_SIMULATOR
