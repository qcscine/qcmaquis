/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/chem/rel/model.h"
#include "dmrg/models/factories/factory.h"

template<class ScalarType>
class RelativisticFactoryClass { };

template<>
class RelativisticFactoryClass<double> {
public:
  static auto getPointerToSimulation(const Lattice& lattice, BaseParameters& parms) {
    typedef std::shared_ptr<model_impl<matrix, U1DG> > impl_ptr;
    throw std::runtime_error("Real-valued relativistic simulation does not make sense!");
    return impl_ptr();
  }
};

template<>
class RelativisticFactoryClass<std::complex<double>> {
public:
  static auto getPointerToSimulation(const Lattice& lattice, BaseParameters& parms) {
    typedef std::shared_ptr<model_impl<cmatrix, U1DG> > impl_ptr;
    return impl_ptr( new rel_qc_model<U1DG>(lattice, parms) );
  }
};

template<class Matrix>
struct coded_model_factory<Matrix, U1DG> {

  static std::shared_ptr<model_impl<Matrix, U1DG> > parse(Lattice const & lattice, BaseParameters & parms)
  {
	  typedef std::shared_ptr<model_impl<Matrix, U1DG> > impl_ptr;
    if (parms["MODEL"] == std::string("relativistic_quantum_chemistry")) {
      if (parms["LATTICE"] != std::string("spinors"))
        throw std::runtime_error("Please use \"LATTICE = spinors\" for relativistic_quantum_chemistry\n");
      return RelativisticFactoryClass<typename Matrix::value_type>::getPointerToSimulation(lattice, parms);
    }
    else {
      throw std::runtime_error("Don't know this model!\n");
      return impl_ptr();
    }
  }
};
