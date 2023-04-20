/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/factories/factory.h"
#include "dmrg/models/vibrational/u1/VibronicModel.hpp"
#include "dmrg/models/vibrational/u1/ExcitonicModel.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    using PointerType = std::shared_ptr<model_impl<Matrix, U1> >;

    /** @brief Factory function for the model class */
    static PointerType parse(Lattice const& lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, U1> > impl_ptr;
        if (parms["MODEL"] == std::string("vibronic")) {
#ifdef DMRG_VIBRONIC
            return impl_ptr( new VibronicModel<Matrix>(lattice, parms));
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else if (parms["MODEL"] == std::string("excitonic")) {
#ifdef DMRG_VIBRONIC
            return impl_ptr( new HolsteinHubbardExcitonicHamiltonian<Matrix>(lattice, parms));
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
