/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
*               2022- Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef FEAST_LAUNCHER
#define FEAST_LAUNCHER

#include <exception>
#include "dmrg/MetaSweepSimulations/FEASTSimulator.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/sim/matrix_types.h"
#include "dmrg/utils/DmrgParameters.h"

/** @brief Exception raised by FEAST NYI functionalities */
struct FEASTException : public std::exception {
  const char* what() const throw () {
    return "FEAST must be run with complex-valued simulations";
  }
};

/**
 * @brief Wrapper class around the construction and running of a FEAST simulation.
 *
 * This class is needed to ensure that a FEAST simulation based on real-valued
 * MPSs and MPOs is not run.
 */
template<class Matrix, class SymmGroup>
class FEASTLauncher {};

template<class SymmGroup>
class FEASTLauncher<matrix, SymmGroup> {
  using ModelType = Model<matrix, SymmGroup>;
  using MPSType = MPS<matrix, SymmGroup>;
  using MPOType = MPO<matrix, SymmGroup>;
  using RetType = std::shared_ptr<std::vector<MPSType>>;
public:
  static RetType runFEASTSimulation(BaseParameters& parms, const ModelType& model,
                                    const Lattice& lattice, const MPOType& mpo)
  {
    throw FEASTException();
  }
};

template<class SymmGroup>
class FEASTLauncher<cmatrix, SymmGroup> {
  using FEASTSimulatorType = FEASTSimulator<SymmGroup>;
  using ModelType = Model<cmatrix, SymmGroup>;
  using MPSType = MPS<cmatrix, SymmGroup>;
  using MPOType = MPO<cmatrix, SymmGroup>;
  using RetType = std::shared_ptr<std::vector<MPSType>>;
public:
  static RetType runFEASTSimulation(BaseParameters& parms, const ModelType& model,
                                    const Lattice& lattice, const MPOType& mpo)
  {
    auto feastSimulator = FEASTSimulatorType(parms, model, lattice, mpo);
    feastSimulator.runFEAST();
    return feastSimulator.getCurrentEigenvalues();
  }
};

#endif // FEAST_LAUNCHER
