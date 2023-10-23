/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
